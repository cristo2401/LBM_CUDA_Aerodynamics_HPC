#include "gpu_lbm_benchmark.hpp"

#include <algorithm>
#include <chrono>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace lbm::gpu_benchmark {
namespace {

constexpr int q_count = 9;

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

__global__ void soa_macro_range_kernel(const float* __restrict__ f,
                                       const unsigned char* __restrict__ solid,
                                       float* __restrict__ rho,
                                       float* __restrict__ ux,
                                       float* __restrict__ uy,
                                       int nx,
                                       int ny,
                                       int y_begin,
                                       int y_end,
                                       float rho0)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = y_begin + blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= nx || y >= y_end || y >= ny) {
        return;
    }

    const int cells = nx * ny;
    const int cell = x + nx * y;
    if (solid[cell]) {
        rho[cell] = rho0;
        ux[cell] = 0.0F;
        uy[cell] = 0.0F;
        return;
    }

    float r = 0.0F;
    float mx = 0.0F;
    float my = 0.0F;
#pragma unroll
    for (int q = 0; q < q_count; ++q) {
        const float fi = f[q * cells + cell];
        r += fi;
        mx += fi * static_cast<float>(d_cx[q]);
        my += fi * static_cast<float>(d_cy[q]);
    }
    if (r <= 1.0e-14F) {
        r = rho0;
        mx = 0.0F;
        my = 0.0F;
    }
    rho[cell] = r;
    ux[cell] = mx / r;
    uy[cell] = my / r;
}

__global__ void soa_collide_stream_shared_kernel(const float* __restrict__ f,
                                                 float* __restrict__ next,
                                                 const float* __restrict__ rho,
                                                 const float* __restrict__ ux,
                                                 const float* __restrict__ uy,
                                                 const unsigned char* __restrict__ solid,
                                                 int nx,
                                                 int ny,
                                                 int y_begin,
                                                 int y_end,
                                                 float omega,
                                                 float rho0)
{
    extern __shared__ float tile[];
    const int tile_w = blockDim.x + 2;
    const int tile_h = blockDim.y + 2;
    const int tile_count = tile_w * tile_h;
    float* s_rho = tile;
    float* s_ux = s_rho + tile_count;
    float* s_uy = s_ux + tile_count;
    float* s_solid = s_uy + tile_count;

    const int base_x = blockIdx.x * blockDim.x;
    const int base_y = y_begin + blockIdx.y * blockDim.y;
    const int cells = nx * ny;

    for (int ly = threadIdx.y; ly < tile_h; ly += blockDim.y) {
        for (int lx = threadIdx.x; lx < tile_w; lx += blockDim.x) {
            const int gx = wrap_coord(base_x + lx - 1, nx);
            const int gy = wrap_coord(base_y + ly - 1, ny);
            const int cell = gx + nx * gy;
            const int local = lx + tile_w * ly;
            s_rho[local] = rho[cell];
            s_ux[local] = ux[cell];
            s_uy[local] = uy[cell];
            s_solid[local] = solid[cell] ? 1.0F : 0.0F;
        }
    }
    __syncthreads();

    const int x = base_x + threadIdx.x;
    const int y = base_y + threadIdx.y;
    if (x >= nx || y >= y_end || y >= ny) {
        return;
    }

    const int cell = x + nx * y;
    const int lx = threadIdx.x + 1;
    const int ly = threadIdx.y + 1;
    const int local = lx + tile_w * ly;

    if (s_solid[local] > 0.5F) {
#pragma unroll
        for (int q = 0; q < q_count; ++q) {
            next[q * cells + cell] = equilibrium_device(q, rho0, 0.0F, 0.0F);
        }
        return;
    }

#pragma unroll
    for (int q = 0; q < q_count; ++q) {
        const int src_x = wrap_coord(x - d_cx[q], nx);
        const int src_y = wrap_coord(y - d_cy[q], ny);
        const int src = src_x + nx * src_y;
        const int slx = lx - d_cx[q];
        const int sly = ly - d_cy[q];
        const int src_local = slx + tile_w * sly;

        int direction = q;
        int read_cell = src;
        float r = s_rho[src_local];
        float vx = s_ux[src_local];
        float vy = s_uy[src_local];

        if (s_solid[src_local] > 0.5F) {
            direction = d_opposite[q];
            read_cell = cell;
            r = s_rho[local];
            vx = s_ux[local];
            vy = s_uy[local];
        }

        const float f_old = f[direction * cells + read_cell];
        const float feq = equilibrium_device(direction, r, vx, vy);
        next[q * cells + cell] = f_old - omega * (f_old - feq);
    }
}

double milliseconds_since(std::chrono::steady_clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start)
        .count();
}

} // namespace

void run_soa_optimized_cuda(const BenchmarkConfig& config,
                            const BenchmarkState& input,
                            BenchmarkState& output,
                            BackendTiming& timing)
{
    output.nx = config.nx;
    output.ny = config.ny;
    const std::size_t cells = input.cells();
    const std::size_t f_bytes = static_cast<std::size_t>(q_count) * cells * sizeof(float);
    const std::size_t scalar_bytes = cells * sizeof(float);
    const std::size_t solid_bytes = cells * sizeof(unsigned char);

    float* d_f = nullptr;
    float* d_next = nullptr;
    float* d_rho = nullptr;
    float* d_ux = nullptr;
    float* d_uy = nullptr;
    unsigned char* d_solid = nullptr;
    cudaStream_t stream0 = nullptr;
    cudaStream_t stream1 = nullptr;

    check_cuda(cudaMalloc(&d_f, f_bytes), "cudaMalloc soa f");
    check_cuda(cudaMalloc(&d_next, f_bytes), "cudaMalloc soa next");
    check_cuda(cudaMalloc(&d_rho, scalar_bytes), "cudaMalloc soa rho");
    check_cuda(cudaMalloc(&d_ux, scalar_bytes), "cudaMalloc soa ux");
    check_cuda(cudaMalloc(&d_uy, scalar_bytes), "cudaMalloc soa uy");
    check_cuda(cudaMalloc(&d_solid, solid_bytes), "cudaMalloc soa solid");
    check_cuda(cudaStreamCreate(&stream0), "create stream0");
    check_cuda(cudaStreamCreate(&stream1), "create stream1");

    try {
        auto start = std::chrono::steady_clock::now();
        check_cuda(cudaMemcpyAsync(d_f, input.f.data(), f_bytes, cudaMemcpyHostToDevice, stream0),
                   "copy soa f H2D");
        check_cuda(cudaMemcpyAsync(d_solid,
                                   input.solid.data(),
                                   solid_bytes,
                                   cudaMemcpyHostToDevice,
                                   stream1),
                   "copy soa solid H2D");
        check_cuda(cudaDeviceSynchronize(), "synchronize soa H2D");
        timing.h2d_ms = milliseconds_since(start);

        start = std::chrono::steady_clock::now();
        const dim3 block(16, 16);
        const dim3 grid((config.nx + block.x - 1) / block.x,
                        (config.ny + block.y - 1) / block.y);
        const std::size_t shared_bytes =
            static_cast<std::size_t>(block.x + 2) *
            static_cast<std::size_t>(block.y + 2) * 4U * sizeof(float);
        const float omega = 1.0F / config.tau;

        for (int step = 0; step < config.steps; ++step) {
            soa_macro_range_kernel<<<grid, block, 0, stream0>>>(
                d_f,
                d_solid,
                d_rho,
                d_ux,
                d_uy,
                config.nx,
                config.ny,
                0,
                config.ny,
                config.rho0);
            soa_collide_stream_shared_kernel<<<grid, block, shared_bytes, stream0>>>(
                d_f,
                d_next,
                d_rho,
                d_ux,
                d_uy,
                d_solid,
                config.nx,
                config.ny,
                0,
                config.ny,
                omega,
                config.rho0);
            std::swap(d_f, d_next);
        }
        soa_macro_range_kernel<<<grid, block, 0, stream0>>>(
            d_f,
            d_solid,
            d_rho,
            d_ux,
            d_uy,
            config.nx,
            config.ny,
            0,
            config.ny,
            config.rho0);
        check_cuda(cudaGetLastError(), "launch soa kernels");
        check_cuda(cudaDeviceSynchronize(), "synchronize soa kernels");
        timing.kernel_ms = milliseconds_since(start);

        output.f.assign(static_cast<std::size_t>(q_count) * cells, 0.0F);
        output.rho.assign(cells, 0.0F);
        output.ux.assign(cells, 0.0F);
        output.uy.assign(cells, 0.0F);
        output.solid.assign(cells, 0);

        start = std::chrono::steady_clock::now();
        check_cuda(cudaMemcpyAsync(output.f.data(), d_f, f_bytes, cudaMemcpyDeviceToHost, stream0),
                   "copy soa f D2H");
        check_cuda(cudaMemcpyAsync(output.rho.data(), d_rho, scalar_bytes, cudaMemcpyDeviceToHost, stream1),
                   "copy soa rho D2H");
        check_cuda(cudaMemcpyAsync(output.ux.data(), d_ux, scalar_bytes, cudaMemcpyDeviceToHost, stream0),
                   "copy soa ux D2H");
        check_cuda(cudaMemcpyAsync(output.uy.data(), d_uy, scalar_bytes, cudaMemcpyDeviceToHost, stream1),
                   "copy soa uy D2H");
        check_cuda(cudaMemcpyAsync(output.solid.data(),
                                   d_solid,
                                   solid_bytes,
                                   cudaMemcpyDeviceToHost,
                                   stream0),
                   "copy soa solid D2H");
        check_cuda(cudaDeviceSynchronize(), "synchronize soa D2H");
        timing.d2h_ms = milliseconds_since(start);
    } catch (...) {
        cudaFree(d_f);
        cudaFree(d_next);
        cudaFree(d_rho);
        cudaFree(d_ux);
        cudaFree(d_uy);
        cudaFree(d_solid);
        cudaStreamDestroy(stream0);
        cudaStreamDestroy(stream1);
        throw;
    }

    check_cuda(cudaFree(d_f), "cudaFree soa f");
    check_cuda(cudaFree(d_next), "cudaFree soa next");
    check_cuda(cudaFree(d_rho), "cudaFree soa rho");
    check_cuda(cudaFree(d_ux), "cudaFree soa ux");
    check_cuda(cudaFree(d_uy), "cudaFree soa uy");
    check_cuda(cudaFree(d_solid), "cudaFree soa solid");
    check_cuda(cudaStreamDestroy(stream0), "destroy stream0");
    check_cuda(cudaStreamDestroy(stream1), "destroy stream1");

    timing.backend = "gpu_soa_optimized";
    timing.layout = "soa_shared_async";
    timing.total_ms = timing.h2d_ms + timing.kernel_ms + timing.d2h_ms;
    timing.mlups = (static_cast<double>(input.cells()) * config.steps) /
        (timing.kernel_ms * 1000.0);
}

} // namespace lbm::gpu_benchmark
