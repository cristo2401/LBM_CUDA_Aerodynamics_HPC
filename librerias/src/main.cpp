#include "case_builder.hpp"
#include "exporter.hpp"
#include "forces.hpp"
#include "lbm_solver.hpp"
#include "observables.hpp"
#include "params.hpp"
#include "utils.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void print_help(const char* program)
{
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Options:\n"
        << "  --nx N                  Domain width (default 420)\n"
        << "  --ny N                  Domain height (default 160)\n"
        << "  --steps N               Time steps (default 2000)\n"
        << "  --frame-stride N        Export one frame every N steps (default 100)\n"
        << "  --force-stride N        Write forces every N steps (default 10)\n"
        << "  --log-stride N          Print progress every N steps (default 100)\n"
        << "  --ramp-steps N          Smooth inlet ramp steps; 0 disables ramp (default 100)\n"
        << "  --u-in U                Inlet velocity in lattice units (default 0.06)\n"
        << "  --tau T                 BGK relaxation time, T > 0.5 (default 0.60)\n"
        << "  --rho0 R                Reference density (default 1.0)\n"
        << "  --obstacle circle|box|airfoil|sdf\n"
        << "                          Obstacle type (default circle)\n"
        << "  --radius R              Circle radius (default 20)\n"
        << "  --hx H                  Box half width (default 30)\n"
        << "  --hy H                  Box half height (default 12)\n"
        << "  --cx X                  Obstacle center x (default 0.25*nx)\n"
        << "  --cy Y                  Obstacle center y (default 0.5*ny)\n"
        << "  --angle-deg A           Box angle, or airfoil AoA in degrees (default 0)\n"
        << "  --chord C               Airfoil chord length (default 80)\n"
        << "  --camber M              Airfoil max camber fraction, e.g. 0.02 (default 0.02)\n"
        << "  --camber-pos P          Airfoil camber position fraction (default 0.40)\n"
        << "  --thickness T           Airfoil max thickness fraction, e.g. 0.12 (default 0.12)\n"
        << "  --airfoil-samples N     Airfoil contour samples per side (default 180)\n"
        << "  --sdf-file PATH         Imported SDF polygon file when --obstacle sdf\n"
        << "  --sdf-scale S           Imported SDF polygon length scale (default 120)\n"
        << "  --boundary open|channel External boundary mode (default open)\n"
        << "  --output DIR            Output directory (default resultados/default)\n"
        << "  --case-name NAME        Case name in metadata\n"
        << "  --no-initial-frame      Do not write frame_000000.bin\n"
        << "  --start-from-rest       Initialize fluid with ux=0 (default)\n"
        << "  --uniform-initial-flow  Initialize fluid with ux=u_in\n"
        << "  --help                  Show this help\n";
}

std::string require_value(int& i, int argc, char** argv, const std::string& option)
{
    if (i + 1 >= argc) {
        throw std::invalid_argument("missing value for " + option);
    }
    return argv[++i];
}

std::size_t parse_size(const std::string& value, const std::string& option)
{
    if (value.empty() || value.front() == '-') {
        throw std::invalid_argument(option + " must be a positive integer");
    }
    std::size_t consumed = 0;
    const auto parsed = std::stoull(value, &consumed);
    if (consumed != value.size()) {
        throw std::invalid_argument(option + " must be a positive integer");
    }
    if (parsed == 0) {
        throw std::invalid_argument(option + " must be greater than zero");
    }
    return static_cast<std::size_t>(parsed);
}

std::size_t parse_nonnegative_size(const std::string& value, const std::string& option)
{
    if (value.empty() || value.front() == '-') {
        throw std::invalid_argument(option + " must be a non-negative integer");
    }
    std::size_t consumed = 0;
    const auto parsed = std::stoull(value, &consumed);
    if (consumed != value.size()) {
        throw std::invalid_argument(option + " must be a non-negative integer");
    }
    return static_cast<std::size_t>(parsed);
}

double parse_double(const std::string& value, const std::string& option)
{
    std::size_t consumed = 0;
    const double parsed = std::stod(value, &consumed);
    if (consumed != value.size()) {
        throw std::invalid_argument(option + " must be a number");
    }
    if (!std::isfinite(parsed)) {
        throw std::invalid_argument(option + " must be finite");
    }
    return parsed;
}

lbm::SimulationParams parse_args(int argc, char** argv, bool& show_help)
{
    lbm::SimulationParams params;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            show_help = true;
            return params;
        }

        if (arg == "--no-initial-frame") {
            params.write_initial_frame = false;
        } else if (arg == "--start-from-rest") {
            params.start_from_rest = true;
        } else if (arg == "--uniform-initial-flow") {
            params.start_from_rest = false;
        } else if (arg == "--nx") {
            params.nx = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--ny") {
            params.ny = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--steps") {
            params.steps = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--frame-stride") {
            params.frame_stride = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--force-stride") {
            params.force_stride = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--log-stride") {
            params.log_stride = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--ramp-steps") {
            params.inlet_ramp_steps = parse_nonnegative_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--u-in") {
            params.u_in = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--tau") {
            params.tau = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--rho0") {
            params.rho0 = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--radius") {
            params.obstacle.radius = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--hx") {
            params.obstacle.half_width = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--hy") {
            params.obstacle.half_height = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--cx") {
            params.obstacle.center_x = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--cy") {
            params.obstacle.center_y = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--angle-deg") {
            params.obstacle.angle_degrees = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--chord") {
            params.obstacle.airfoil_chord = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--camber") {
            params.obstacle.airfoil_max_camber = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--camber-pos") {
            params.obstacle.airfoil_camber_position = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--thickness") {
            params.obstacle.airfoil_max_thickness = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--airfoil-samples") {
            params.obstacle.airfoil_samples = parse_size(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--sdf-file") {
            params.obstacle.sdf_file = require_value(i, argc, argv, arg);
        } else if (arg == "--sdf-scale") {
            params.obstacle.sdf_scale = parse_double(require_value(i, argc, argv, arg), arg);
        } else if (arg == "--output") {
            params.output_dir = require_value(i, argc, argv, arg);
        } else if (arg == "--case-name") {
            params.case_name = require_value(i, argc, argv, arg);
        } else if (arg == "--obstacle") {
            const std::string value = require_value(i, argc, argv, arg);
            if (value == "circle") {
                params.obstacle.type = lbm::ObstacleType::Circle;
            } else if (value == "box") {
                params.obstacle.type = lbm::ObstacleType::Box;
            } else if (value == "airfoil") {
                params.obstacle.type = lbm::ObstacleType::Airfoil;
            } else if (value == "sdf") {
                params.obstacle.type = lbm::ObstacleType::ImportedSdf;
            } else {
                throw std::invalid_argument("--obstacle must be circle, box, airfoil, or sdf");
            }
        } else if (arg == "--boundary") {
            const std::string value = require_value(i, argc, argv, arg);
            if (value == "open") {
                params.boundary_mode = lbm::BoundaryMode::OpenFlow;
            } else if (value == "channel") {
                params.boundary_mode = lbm::BoundaryMode::Channel;
            } else {
                throw std::invalid_argument("--boundary must be open or channel");
            }
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }

    return params;
}

void print_summary(const lbm::SimulationParams& params)
{
    std::cout << std::fixed << std::setprecision(4)
              << "case=" << params.case_name
              << " grid=" << params.nx << 'x' << params.ny
              << " steps=" << params.steps
              << " u_in=" << params.u_in
              << " tau=" << params.tau
              << " nu=" << params.viscosity()
              << " Re=" << params.reynolds()
              << " init=" << (params.start_from_rest ? "rest" : "uniform")
              << " ramp_steps=" << params.inlet_ramp_steps
              << " output=" << params.output_dir
              << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        bool show_help = false;
        lbm::SimulationParams params = parse_args(argc, argv, show_help);
        if (show_help) {
            print_help(argv[0]);
            return EXIT_SUCCESS;
        }

        params.validate();
        print_summary(params);

        lbm::case_builder::Case sim_case = lbm::case_builder::build(params);
        std::cout << "cells: fluid=" << sim_case.fluid_cells
                  << " solid=" << sim_case.solid_cells << '\n';
        for (const std::string& warning : sim_case.warnings) {
            std::cerr << "warning: " << warning << '\n';
        }

        lbm::LBMSolver solver(params, std::move(sim_case.field));
        solver.initialize();

        lbm::Exporter exporter(params);
        exporter.prepare(solver.field());

        if (params.write_initial_frame) {
            exporter.write_frame(solver.field(), 0);
        }
        exporter.append_force(lbm::forces::compute(solver.field(), params, 0));
        exporter.append_diagnostics(lbm::observables::compute(solver.field(), params, 0));

        lbm::utils::Timer timer;
        lbm::forces::ForceSample last_force;
        lbm::observables::Sample last_diag;
        for (std::size_t step = 1; step <= params.steps; ++step) {
            solver.step();

            bool force_valid = false;
            bool diag_valid = false;
            if (step % params.force_stride == 0) {
                last_force = lbm::forces::compute(solver.field(), params, step);
                last_diag = lbm::observables::compute(solver.field(), params, step);
                exporter.append_force(last_force);
                exporter.append_diagnostics(last_diag);
                force_valid = true;
                diag_valid = true;
            }

            if (step % params.frame_stride == 0) {
                exporter.write_frame(solver.field(), step);
            }

            if (step % params.log_stride == 0 || step == params.steps) {
                if (!force_valid) {
                    last_force = lbm::forces::compute(solver.field(), params, step);
                }
                if (!diag_valid) {
                    last_diag = lbm::observables::compute(solver.field(), params, step);
                }
                const double elapsed = timer.elapsed_seconds();
                const double mlups = elapsed > 0.0
                    ? (static_cast<double>(params.nx * params.ny * step) / elapsed) / 1.0e6
                    : 0.0;

                std::cout << "step=" << step << '/' << params.steps
                          << " Cd=" << last_force.cd
                          << " Cl=" << last_force.cl
                          << " maxU=" << last_diag.max_speed
                          << " rho=[" << last_diag.min_rho << ',' << last_diag.max_rho << ']'
                          << " frames=" << exporter.frame_count()
                          << " MLUPS=" << mlups
                          << '\n';

                if (!last_diag.finite) {
                    throw std::runtime_error("non-finite value detected in simulation fields.");
                }
            }
        }

        std::cout << "done. metadata: " << (params.output_dir + "/metadata.json") << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return EXIT_FAILURE;
    }
}
