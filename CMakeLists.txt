cmake_minimum_required(VERSION 3.16)

project(lbm_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_OPENMP "Enable OpenMP parallel loops when available" ON)

add_executable(lbm2d
    librerias/src/main.cpp
    librerias/src/case_builder.cpp
    librerias/src/lbm_solver.cpp
    librerias/src/lbm_utils.cpp
    librerias/src/sdf.cpp
    librerias/src/geometry.cpp
    librerias/src/boundary.cpp
    librerias/src/forces.cpp
    librerias/src/observables.cpp
    librerias/src/exporter.cpp
    librerias/src/utils.cpp
)

target_include_directories(lbm2d PRIVATE librerias/include)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(lbm2d PRIVATE -Wall -Wextra -Wpedantic -O3)
endif()

if(ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(lbm2d PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(lbm2d PRIVATE LBM_USE_OPENMP)
    endif()
endif()

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    add_executable(lbm_gpu_benchmark
        librerias/benchmarks/gpu_lbm_benchmark.cpp
        librerias/cuda/lbm_aos_naive.cu
        librerias/cuda/lbm_soa_optimized.cu
    )
    target_include_directories(lbm_gpu_benchmark PRIVATE librerias/include)
    set_target_properties(lbm_gpu_benchmark PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF
    )
endif()
