#include "gpu_lbm_benchmark.hpp"

#include <chrono>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace lbm::gpu_benchmark {
namespace {

constexpr int q_count = 9;

struct AoSCell {
    float f[q_count];
    float rho;
    float ux;
    float uy;
    unsigned char solid;
};

__constant__ int d_cx[q_count] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
__constant__ int d_cy[q_count] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
__constant__ int d_opposite[q_count] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
__constant__ float d_w[q_count] = {
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

void check_cuda(cudaError_t status, const char* what)
{
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
    }
}

__device__ __forceinline__ int wrap_coord(int value, int limit)
{
    if (value < 0) {
        return value + limit;
    }
    if (value >= limit) {
        return value - limit;
    }
    return value;
}

__device__ __forceinline__ float equilibrium_device(int direction,
                                                    float rho,
                                                    float ux,
                                                    float uy)
{
    const float cu = static_cast<float>(d_cx[direction]) * ux +
                     static_cast<float>(d_cy[direction]) * uy;
    const float u2 = ux * ux + uy * uy;
    return d_w[direction] * rho *
           (1.0F + 3.0F * cu + 4.5F * cu * cu - 1.5F * u2);
}

__global__ void aos_macro_kernel(AoSCell* cells, int total_cells, float rho0)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= total_cells) {
        return;
    }

    AoSCell& c = cells[cell];
    if (c.solid) {
        c.rho = rho0;
        c.ux = 0.0F;
        c.uy = 0.0F;
        return;
    }

    float rho = 0.0F;
    float mx = 0.0F;
    float my = 0.0F;
    for (int q = 0; q < q_count; ++q) {
        const float fi = c.f[q];
        rho += fi;
        mx += fi * static_cast<float>(d_cx[q]);
        my += fi * static_cast<float>(d_cy[q]);
    }
    if (rho <= 1.0e-14F) {
        rho = rho0;
        mx = 0.0F;
        my = 0.0F;
    }
    c.rho = rho;
    c.ux = mx / rho;
    c.uy = my / rho;
}

__device__ __forceinline__ float collide_aos(const AoSCell* cells,
                                             int direction,
                                             int cell,
                                             float omega)
{
    const AoSCell& c = cells[cell];
    const float f_old = c.f[direction];
    const float feq = equilibrium_device(direction, c.rho, c.ux, c.uy);
    return f_old - omega * (f_old - feq);
}

__global__ void aos_collide_stream_kernel(const AoSCell* current,
                                          AoSCell* next,
                                          int nx,
                                          int ny,
                                          float omega,
                                          float rho0)
{
    const int cell = blockIdx.x * blockDim.x + threadIdx.x;
    const int total_cells = nx * ny;
    if (cell >= total_cells) {
        return;
    }

    const int x = cell % nx;
    const int y = cell / nx;
    const AoSCell& c = current[cell];
    AoSCell out = {};
    out.solid = c.solid;

    if (c.solid) {
        out.rho = rho0;
        out.ux = 0.0F;
        out.uy = 0.0F;
        for (int q = 0; q < q_count; ++q) {
            out.f[q] = equilibrium_device(q, rho0, 0.0F, 0.0F);
        }
        next[cell] = out;
        return;
    }

    for (int q = 0; q < q_count; ++q) {
        const int sx = wrap_coord(x - d_cx[q], nx);
        const int sy = wrap_coord(y - d_cy[q], ny);
        const int src = sx + nx * sy;
        out.f[q] = current[src].solid
            ? collide_aos(current, d_opposite[q], cell, omega)
            : collide_aos(current, q, src, omega);
    }

    out.rho = c.rho;
    out.ux = c.ux;
    out.uy = c.uy;
    next[cell] = out;
}

void pack_aos(const BenchmarkState& input, std::vector<AoSCell>& cells)
{
    const std::size_t count = input.cells();
    cells.assign(count, {});
    for (std::size_t cell = 0; cell < count; ++cell) {
        for (int q = 0; q < q_count; ++q) {
            cells[cell].f[q] = input.f[static_cast<std::size_t>(q) * count + cell];
        }
        cells[cell].rho = input.rho[cell];
        cells[cell].ux = input.ux[cell];
        cells[cell].uy = input.uy[cell];
        cells[cell].solid = input.solid[cell];
    }
}

void unpack_aos(const std::vector<AoSCell>& cells, BenchmarkState& output)
{
    const std::size_t count = output.cells();
    output.f.assign(static_cast<std::size_t>(q_count) * count, 0.0F);
    output.rho.assign(count, 0.0F);
    output.ux.assign(count, 0.0F);
    output.uy.assign(count, 0.0F);
    output.solid.assign(count, 0);
    for (std::size_t cell = 0; cell < count; ++cell) {
        for (int q = 0; q < q_count; ++q) {
            output.f[static_cast<std::size_t>(q) * count + cell] = cells[cell].f[q];
        }
        output.rho[cell] = cells[cell].rho;
        output.ux[cell] = cells[cell].ux;
        output.uy[cell] = cells[cell].uy;
        output.solid[cell] = cells[cell].solid;
    }
}

double milliseconds_since(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start)
        .count();
}

} // namespace

void run_aos_naive_cuda(const BenchmarkConfig& config,
                        const BenchmarkState& input,
                        BenchmarkState& output,
                        BackendTiming& timing)
{
    output.nx = config.nx;
    output.ny = config.ny;

    const int cells = config.nx * config.ny;
    const std::size_t bytes = static_cast<std::size_t>(cells) * sizeof(AoSCell);
    std::vector<AoSCell> host_cells;
    std::vector<AoSCell> host_result(static_cast<std::size_t>(cells));
    pack_aos(input, host_cells);

    AoSCell* d_current = nullptr;
    AoSCell* d_next = nullptr;
    check_cuda(cudaMalloc(&d_current, bytes), "cudaMalloc aos current");
    check_cuda(cudaMalloc(&d_next, bytes), "cudaMalloc aos next");

    try {
        auto start = std::chrono::steady_clock::now();
        check_cuda(cudaMemcpy(d_current, host_cells.data(), bytes, cudaMemcpyHostToDevice),
                   "cudaMemcpy aos H2D");
        check_cuda(cudaDeviceSynchronize(), "cuda synchronize aos H2D");
        timing.h2d_ms = milliseconds_since(start);

        const int threads = 128;
        const int blocks = (cells + threads - 1) / threads;
        const float omega = 1.0F / config.tau;

        start = std::chrono::steady_clock::now();
        for (int step = 0; step < config.steps; ++step) {
            aos_macro_kernel<<<blocks, threads>>>(d_current, cells, config.rho0);
            aos_collide_stream_kernel<<<blocks, threads>>>(
                d_current, d_next, config.nx, config.ny, omega, config.rho0);
            std::swap(d_current, d_next);
        }
        aos_macro_kernel<<<blocks, threads>>>(d_current, cells, config.rho0);
        check_cuda(cudaGetLastError(), "launch aos kernels");
        check_cuda(cudaDeviceSynchronize(), "cuda synchronize aos kernels");
        timing.kernel_ms = milliseconds_since(start);

        start = std::chrono::steady_clock::now();
        check_cuda(cudaMemcpy(host_result.data(), d_current, bytes, cudaMemcpyDeviceToHost),
                   "cudaMemcpy aos D2H");
        check_cuda(cudaDeviceSynchronize(), "cuda synchronize aos D2H");
        timing.d2h_ms = milliseconds_since(start);
    } catch (...) {
        cudaFree(d_current);
        cudaFree(d_next);
        throw;
    }

    check_cuda(cudaFree(d_current), "cudaFree aos current");
    check_cuda(cudaFree(d_next), "cudaFree aos next");

    unpack_aos(host_result, output);
    timing.backend = "gpu_aos_naive";
    timing.layout = "aos";
    timing.total_ms = timing.h2d_ms + timing.kernel_ms + timing.d2h_ms;
    timing.mlups = (static_cast<double>(input.cells()) * config.steps) /
        (timing.kernel_ms * 1000.0);
}

} // namespace lbm::gpu_benchmark
