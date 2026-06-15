#include "gpu_lbm_benchmark.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace lbm::gpu_benchmark {
namespace {

constexpr int q_count = 9;
constexpr int cx[q_count] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr int cy[q_count] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr int opposite[q_count] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
constexpr float w[q_count] = {
    4.0F / 9.0F,
    1.0F / 9.0F,
    1.0F / 9.0F,
    1.0F / 9.0F,
    1.0F / 9.0F,
    1.0F / 36.0F,
    1.0F / 36.0F,
    1.0F / 36.0F,
    1.0F / 36.0F
};

struct Difference {
    double max_abs_rho = 0.0;
    double mean_abs_rho = 0.0;
    double max_abs_ux = 0.0;
    double mean_abs_ux = 0.0;
    double max_abs_uy = 0.0;
    double mean_abs_uy = 0.0;
    double max_abs_f = 0.0;
    double mean_abs_f = 0.0;
};

float equilibrium(int direction, float rho, float ux, float uy)
{
    const float cu = static_cast<float>(cx[direction]) * ux +
                     static_cast<float>(cy[direction]) * uy;
    const float u2 = ux * ux + uy * uy;
    return w[direction] * rho *
           (1.0F + 3.0F * cu + 4.5F * cu * cu - 1.5F * u2);
}

std::size_t index(int x, int y, int nx)
{
    return static_cast<std::size_t>(x) + static_cast<std::size_t>(nx) *
        static_cast<std::size_t>(y);
}

std::size_t dist_index(int direction, std::size_t cell, std::size_t cells)
{
    return static_cast<std::size_t>(direction) * cells + cell;
}

BenchmarkState make_initial_state(const BenchmarkConfig& config)
{
    BenchmarkState state;
    state.nx = config.nx;
    state.ny = config.ny;
    const std::size_t cells = state.cells();
    state.f.assign(static_cast<std::size_t>(q_count) * cells, 0.0F);
    state.rho.assign(cells, config.rho0);
    state.ux.assign(cells, 0.0F);
    state.uy.assign(cells, 0.0F);
    state.solid.assign(cells, 0);

    if (!config.solid_frame.empty()) {
        std::ifstream in(config.solid_frame, std::ios::binary);
        if (!in) {
            throw std::runtime_error("could not open solid frame: " + config.solid_frame);
        }
        const std::streamoff offset =
            static_cast<std::streamoff>(config.solid_field_index) *
            static_cast<std::streamoff>(cells) *
            static_cast<std::streamoff>(sizeof(float));
        in.seekg(offset);
        if (!in) {
            throw std::runtime_error("could not seek solid field in: " + config.solid_frame);
        }
        std::vector<float> solid_values(cells, 0.0F);
        in.read(reinterpret_cast<char*>(solid_values.data()),
                static_cast<std::streamsize>(solid_values.size() * sizeof(float)));
        if (!in) {
            throw std::runtime_error("could not read solid field from: " + config.solid_frame);
        }
        for (std::size_t cell = 0; cell < cells; ++cell) {
            state.solid[cell] = solid_values[cell] > 0.5F ? 1U : 0U;
        }
    } else {
        const float cx0 = 0.40F * static_cast<float>(config.nx);
        const float cy0 = 0.50F * static_cast<float>(config.ny);
        const float r2 = config.obstacle_radius * config.obstacle_radius;
        for (int y = 0; y < config.ny; ++y) {
            for (int x = 0; x < config.nx; ++x) {
                const std::size_t cell = index(x, y, config.nx);
                const float dx = static_cast<float>(x) - cx0;
                const float dy = static_cast<float>(y) - cy0;
                state.solid[cell] = dx * dx + dy * dy <= r2 ? 1U : 0U;
            }
        }
    }

    const float pi = std::acos(-1.0F);

    for (int y = 0; y < config.ny; ++y) {
        for (int x = 0; x < config.nx; ++x) {
            const std::size_t cell = index(x, y, config.nx);
            const bool solid = state.solid[cell] != 0;

            const float wave = std::sin(2.0F * pi * static_cast<float>(x) /
                                        static_cast<float>(config.nx));
            const float shear = std::cos(2.0F * pi * static_cast<float>(y) /
                                         static_cast<float>(config.ny));
            const float ux = solid ? 0.0F : config.u_in * (1.0F + 0.02F * wave);
            const float uy = solid ? 0.0F : 0.004F * config.u_in * shear;

            state.rho[cell] = config.rho0;
            state.ux[cell] = ux;
            state.uy[cell] = uy;
            for (int q = 0; q < q_count; ++q) {
                state.f[dist_index(q, cell, cells)] = equilibrium(q, config.rho0, ux, uy);
            }
        }
    }

    return state;
}

void compute_macroscopic(BenchmarkState& state, float rho0)
{
    const std::size_t cells = state.cells();
    for (std::size_t cell = 0; cell < cells; ++cell) {
        if (state.solid[cell]) {
            state.rho[cell] = rho0;
            state.ux[cell] = 0.0F;
            state.uy[cell] = 0.0F;
            continue;
        }

        float rho = 0.0F;
        float mx = 0.0F;
        float my = 0.0F;
        for (int q = 0; q < q_count; ++q) {
            const float fi = state.f[dist_index(q, cell, cells)];
            rho += fi;
            mx += fi * static_cast<float>(cx[q]);
            my += fi * static_cast<float>(cy[q]);
        }
        if (rho <= 1.0e-14F) {
            rho = rho0;
            mx = 0.0F;
            my = 0.0F;
        }
        state.rho[cell] = rho;
        state.ux[cell] = mx / rho;
        state.uy[cell] = my / rho;
    }
}

float collide_value(const BenchmarkState& state,
                    int direction,
                    std::size_t cell,
                    float omega)
{
    const std::size_t cells = state.cells();
    const float f_old = state.f[dist_index(direction, cell, cells)];
    const float feq = equilibrium(direction, state.rho[cell], state.ux[cell], state.uy[cell]);
    return f_old - omega * (f_old - feq);
}

void cpu_step(BenchmarkState& state, std::vector<float>& next, float omega, float rho0)
{
    const int nx = state.nx;
    const int ny = state.ny;
    const std::size_t cells = state.cells();
    std::fill(next.begin(), next.end(), 0.0F);

    compute_macroscopic(state, rho0);

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const std::size_t cell = index(x, y, nx);
            if (state.solid[cell]) {
                for (int q = 0; q < q_count; ++q) {
                    next[dist_index(q, cell, cells)] = equilibrium(q, rho0, 0.0F, 0.0F);
                }
                continue;
            }

            for (int q = 0; q < q_count; ++q) {
                const int sx = (x - cx[q] + nx) % nx;
                const int sy = (y - cy[q] + ny) % ny;
                const std::size_t src = index(sx, sy, nx);
                const float streamed = state.solid[src]
                    ? collide_value(state, opposite[q], cell, omega)
                    : collide_value(state, q, src, omega);
                next[dist_index(q, cell, cells)] = streamed;
            }
        }
    }

    state.f.swap(next);
}

void run_cpu_reference(const BenchmarkConfig& config,
                       const BenchmarkState& input,
                       BenchmarkState& output,
                       BackendTiming& timing)
{
    output = input;
    std::vector<float> next(output.f.size(), 0.0F);
    const float omega = 1.0F / config.tau;

    const auto start = std::chrono::steady_clock::now();
    for (int step = 0; step < config.steps; ++step) {
        cpu_step(output, next, omega, config.rho0);
    }
    compute_macroscopic(output, config.rho0);
    const auto stop = std::chrono::steady_clock::now();

    timing.backend = "cpu_reference";
    timing.layout = "soa_scalar";
    timing.kernel_ms = std::chrono::duration<double, std::milli>(stop - start).count();
    timing.total_ms = timing.kernel_ms;
    timing.mlups = (static_cast<double>(output.cells()) * config.steps) /
        (timing.kernel_ms * 1000.0);
}

Difference compare_states(const BenchmarkState& reference, const BenchmarkState& candidate)
{
    if (reference.nx != candidate.nx || reference.ny != candidate.ny ||
        reference.f.size() != candidate.f.size()) {
        throw std::runtime_error("cannot compare states with different shapes.");
    }

    Difference diff;
    const std::size_t cells = reference.cells();
    for (std::size_t cell = 0; cell < cells; ++cell) {
        const double drho = std::abs(static_cast<double>(reference.rho[cell]) -
                                     static_cast<double>(candidate.rho[cell]));
        const double dux = std::abs(static_cast<double>(reference.ux[cell]) -
                                    static_cast<double>(candidate.ux[cell]));
        const double duy = std::abs(static_cast<double>(reference.uy[cell]) -
                                    static_cast<double>(candidate.uy[cell]));
        diff.max_abs_rho = std::max(diff.max_abs_rho, drho);
        diff.max_abs_ux = std::max(diff.max_abs_ux, dux);
        diff.max_abs_uy = std::max(diff.max_abs_uy, duy);
        diff.mean_abs_rho += drho;
        diff.mean_abs_ux += dux;
        diff.mean_abs_uy += duy;
    }

    for (std::size_t i = 0; i < reference.f.size(); ++i) {
        const double df = std::abs(static_cast<double>(reference.f[i]) -
                                   static_cast<double>(candidate.f[i]));
        diff.max_abs_f = std::max(diff.max_abs_f, df);
        diff.mean_abs_f += df;
    }

    diff.mean_abs_rho /= static_cast<double>(cells);
    diff.mean_abs_ux /= static_cast<double>(cells);
    diff.mean_abs_uy /= static_cast<double>(cells);
    diff.mean_abs_f /= static_cast<double>(reference.f.size());
    return diff;
}

void write_csv_row(std::ofstream& out,
                   const BenchmarkConfig& config,
                   const BackendTiming& timing,
                   const Difference& diff)
{
    out << timing.backend << ','
        << timing.layout << ','
        << config.nx << ','
        << config.ny << ','
        << config.steps << ','
        << static_cast<std::size_t>(config.nx) * static_cast<std::size_t>(config.ny) << ','
        << timing.h2d_ms << ','
        << timing.kernel_ms << ','
        << timing.d2h_ms << ','
        << timing.total_ms << ','
        << timing.mlups << ','
        << diff.max_abs_rho << ','
        << diff.mean_abs_rho << ','
        << diff.max_abs_ux << ','
        << diff.mean_abs_ux << ','
        << diff.max_abs_uy << ','
        << diff.mean_abs_uy << ','
        << diff.max_abs_f << ','
        << diff.mean_abs_f << '\n';
}

void write_outputs(const BenchmarkConfig& config,
                   const BackendTiming& cpu_timing,
                   const BackendTiming& aos_timing,
                   const BackendTiming& soa_timing,
                   const Difference& aos_diff,
                   const Difference& soa_diff)
{
    std::filesystem::create_directories(config.output_dir);
    const std::filesystem::path out_dir(config.output_dir);
    const auto csv_path = out_dir / "gpu_kernel_timings.csv";
    const auto json_path = out_dir / "gpu_kernel_comparison.json";

    std::ofstream csv(csv_path, std::ios::trunc);
    csv << std::setprecision(9);
    csv << "backend,layout,nx,ny,steps,cells,h2d_ms,kernel_ms,d2h_ms,total_ms,mlups,"
           "max_abs_rho,mean_abs_rho,max_abs_ux,mean_abs_ux,max_abs_uy,mean_abs_uy,"
           "max_abs_f,mean_abs_f\n";
    write_csv_row(csv, config, cpu_timing, {});
    write_csv_row(csv, config, aos_timing, aos_diff);
    write_csv_row(csv, config, soa_timing, soa_diff);

    std::ofstream json(json_path, std::ios::trunc);
    json << std::setprecision(9);
    json << "{\n";
    json << "  \"nx\": " << config.nx << ",\n";
    json << "  \"ny\": " << config.ny << ",\n";
    json << "  \"steps\": " << config.steps << ",\n";
    json << "  \"rho0\": " << config.rho0 << ",\n";
    json << "  \"u_in\": " << config.u_in << ",\n";
    json << "  \"tau\": " << config.tau << ",\n";
    json << "  \"solid_frame\": \"" << config.solid_frame << "\",\n";
    json << "  \"solid_field_index\": " << config.solid_field_index << ",\n";
    json << "  \"timings_csv\": \"" << csv_path.string() << "\",\n";
    json << "  \"backends\": [\n";
    const BackendTiming timings[3] = {cpu_timing, aos_timing, soa_timing};
    for (int i = 0; i < 3; ++i) {
        json << "    {\"backend\": \"" << timings[i].backend
             << "\", \"layout\": \"" << timings[i].layout
             << "\", \"h2d_ms\": " << timings[i].h2d_ms
             << ", \"kernel_ms\": " << timings[i].kernel_ms
             << ", \"d2h_ms\": " << timings[i].d2h_ms
             << ", \"total_ms\": " << timings[i].total_ms
             << ", \"mlups\": " << timings[i].mlups << "}";
        json << (i == 2 ? "\n" : ",\n");
    }
    json << "  ],\n";
    json << "  \"comparison_vs_cpu\": {\n";
    json << "    \"gpu_aos_naive\": {\"max_abs_rho\": " << aos_diff.max_abs_rho
         << ", \"mean_abs_rho\": " << aos_diff.mean_abs_rho
         << ", \"max_abs_ux\": " << aos_diff.max_abs_ux
         << ", \"mean_abs_ux\": " << aos_diff.mean_abs_ux
         << ", \"max_abs_uy\": " << aos_diff.max_abs_uy
         << ", \"mean_abs_uy\": " << aos_diff.mean_abs_uy
         << ", \"max_abs_f\": " << aos_diff.max_abs_f
         << ", \"mean_abs_f\": " << aos_diff.mean_abs_f << "},\n";
    json << "    \"gpu_soa_optimized\": {\"max_abs_rho\": " << soa_diff.max_abs_rho
         << ", \"mean_abs_rho\": " << soa_diff.mean_abs_rho
         << ", \"max_abs_ux\": " << soa_diff.max_abs_ux
         << ", \"mean_abs_ux\": " << soa_diff.mean_abs_ux
         << ", \"max_abs_uy\": " << soa_diff.max_abs_uy
         << ", \"mean_abs_uy\": " << soa_diff.mean_abs_uy
         << ", \"max_abs_f\": " << soa_diff.max_abs_f
         << ", \"mean_abs_f\": " << soa_diff.mean_abs_f << "}\n";
    json << "  }\n";
    json << "}\n";
}

BenchmarkConfig parse_args(int argc, char** argv)
{
    BenchmarkConfig config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) {
                throw std::invalid_argument("missing value for " + option);
            }
            return argv[++i];
        };

        if (arg == "--nx") {
            config.nx = std::stoi(require_value(arg));
        } else if (arg == "--ny") {
            config.ny = std::stoi(require_value(arg));
        } else if (arg == "--steps") {
            config.steps = std::stoi(require_value(arg));
        } else if (arg == "--u-in") {
            config.u_in = std::stof(require_value(arg));
        } else if (arg == "--tau") {
            config.tau = std::stof(require_value(arg));
        } else if (arg == "--radius") {
            config.obstacle_radius = std::stof(require_value(arg));
        } else if (arg == "--output") {
            config.output_dir = require_value(arg);
        } else if (arg == "--solid-frame") {
            config.solid_frame = require_value(arg);
        } else if (arg == "--solid-field-index") {
            config.solid_field_index = std::stoi(require_value(arg));
        } else if (arg == "--help" || arg == "-h") {
            std::cout
                << "Usage: lbm_gpu_benchmark [--nx N] [--ny N] [--steps N]\n"
                << "                         [--u-in U] [--tau T] [--radius R]\n"
                << "                         [--solid-frame frame.bin]\n"
                << "                         [--solid-field-index I]\n"
                << "                         [--output DIR]\n";
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    if (config.nx < 16 || config.ny < 16) {
        throw std::invalid_argument("nx and ny must be at least 16.");
    }
    if (config.steps <= 0) {
        throw std::invalid_argument("steps must be positive.");
    }
    if (config.tau <= 0.5F) {
        throw std::invalid_argument("tau must be greater than 0.5.");
    }
    if (config.solid_field_index < 0) {
        throw std::invalid_argument("solid field index must be non-negative.");
    }
    return config;
}

void print_row(const BackendTiming& timing)
{
    std::cout << std::left << std::setw(22) << timing.backend
              << std::setw(14) << timing.layout
              << " kernel_ms=" << std::setw(10) << std::right << timing.kernel_ms
              << " total_ms=" << std::setw(10) << timing.total_ms
              << " MLUPS=" << timing.mlups << '\n';
}

} // namespace
} // namespace lbm::gpu_benchmark

int main(int argc, char** argv)
{
    try {
        using namespace lbm::gpu_benchmark;

        const BenchmarkConfig config = parse_args(argc, argv);
        const BenchmarkState initial = make_initial_state(config);

        BenchmarkState cpu_state;
        BenchmarkState aos_state;
        BenchmarkState soa_state;
        BackendTiming cpu_timing;
        BackendTiming aos_timing;
        BackendTiming soa_timing;

        run_cpu_reference(config, initial, cpu_state, cpu_timing);
        run_aos_naive_cuda(config, initial, aos_state, aos_timing);
        run_soa_optimized_cuda(config, initial, soa_state, soa_timing);

        const Difference aos_diff = compare_states(cpu_state, aos_state);
        const Difference soa_diff = compare_states(cpu_state, soa_state);
        write_outputs(config, cpu_timing, aos_timing, soa_timing, aos_diff, soa_diff);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "LBM kernel benchmark " << config.nx << 'x' << config.ny
                  << " steps=" << config.steps << '\n';
        print_row(cpu_timing);
        print_row(aos_timing);
        print_row(soa_timing);
        std::cout << "AOS vs CPU max_abs_rho=" << aos_diff.max_abs_rho
                  << " max_abs_ux=" << aos_diff.max_abs_ux
                  << " max_abs_uy=" << aos_diff.max_abs_uy << '\n';
        std::cout << "SOA vs CPU max_abs_rho=" << soa_diff.max_abs_rho
                  << " max_abs_ux=" << soa_diff.max_abs_ux
                  << " max_abs_uy=" << soa_diff.max_abs_uy << '\n';
        std::cout << "wrote " << config.output_dir << "/gpu_kernel_timings.csv and "
                  << config.output_dir << "/gpu_kernel_comparison.json\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
