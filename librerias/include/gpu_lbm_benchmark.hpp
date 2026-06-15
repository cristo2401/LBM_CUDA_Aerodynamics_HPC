#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace lbm::gpu_benchmark {

struct BenchmarkConfig {
    int nx = 512;
    int ny = 192;
    int steps = 300;
    float rho0 = 1.0F;
    float u_in = 0.04F;
    float tau = 0.60F;
    float obstacle_radius = 32.0F;
    std::string output_dir = "data";
    std::string solid_frame;
    int solid_field_index = 5;
};

struct BenchmarkState {
    int nx = 0;
    int ny = 0;
    std::vector<float> f;
    std::vector<float> rho;
    std::vector<float> ux;
    std::vector<float> uy;
    std::vector<unsigned char> solid;

    [[nodiscard]] std::size_t cells() const
    {
        return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    }
};

struct BackendTiming {
    std::string backend;
    std::string layout;
    double h2d_ms = 0.0;
    double kernel_ms = 0.0;
    double d2h_ms = 0.0;
    double total_ms = 0.0;
    double mlups = 0.0;
};

void run_aos_naive_cuda(const BenchmarkConfig& config,
                        const BenchmarkState& input,
                        BenchmarkState& output,
                        BackendTiming& timing);

void run_soa_optimized_cuda(const BenchmarkConfig& config,
                            const BenchmarkState& input,
                            BenchmarkState& output,
                            BackendTiming& timing);

} // namespace lbm::gpu_benchmark
