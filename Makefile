CXX ?= g++
NVCC ?= nvcc
BUILD_DIR ?= ejecuciones/build
SRC_DIR := librerias/src
INC_DIR := librerias/include
CUDA_DIR := librerias/cuda
BENCH_DIR := librerias/benchmarks
GEOM_DIR := librerias/geometries
RESULTS_DIR := resultados
TARGET ?= $(BUILD_DIR)/lbm2d
CUDA_TARGET ?= $(BUILD_DIR)/lbm_gpu_benchmark
OPENMP ?= 1

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/obj/%.o,$(SRC))
DEP := $(OBJ:.o=.d)

CPPFLAGS += -I$(INC_DIR)
CXXFLAGS += -std=c++17 -O3 -Wall -Wextra -Wpedantic
CUDAFLAGS += -std=c++17 -O3 -I$(INC_DIR)
DEPFLAGS += -MMD -MP

ifeq ($(OPENMP),1)
CXXFLAGS += -fopenmp -DLBM_USE_OPENMP
LDLIBS += -fopenmp
endif

.PHONY: all run smoke vehicle-smoke airfoil sports-car formula1 gpu-benchmark gpu-benchmark-f1 clean

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $(OBJ) $(LDLIBS) -o $@

$(BUILD_DIR)/obj/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(CUDA_TARGET): $(BENCH_DIR)/gpu_lbm_benchmark.cpp $(CUDA_DIR)/lbm_aos_naive.cu $(CUDA_DIR)/lbm_soa_optimized.cu $(INC_DIR)/gpu_lbm_benchmark.hpp
	@mkdir -p $(dir $@)
	$(NVCC) $(CUDAFLAGS) $(BENCH_DIR)/gpu_lbm_benchmark.cpp $(CUDA_DIR)/lbm_aos_naive.cu $(CUDA_DIR)/lbm_soa_optimized.cu -o $@

run: $(TARGET)
	./$(TARGET)

smoke: $(TARGET)
	./$(TARGET) --nx 96 --ny 48 --steps 50 --frame-stride 25 --force-stride 5 --log-stride 25 --radius 6 --output $(RESULTS_DIR)/smoke_test --case-name smoke_test

vehicle-smoke: $(TARGET)
	./$(TARGET) --nx 180 --ny 80 --steps 30 --frame-stride 10 --force-stride 5 --log-stride 15 --ramp-steps 10 --obstacle sdf --sdf-file $(GEOM_DIR)/sports_car.sdfpoly --sdf-scale 55 --u-in 0.025 --tau 0.62 --output $(RESULTS_DIR)/smoke_sports_car --case-name smoke_sports_car
	./$(TARGET) --nx 200 --ny 80 --steps 30 --frame-stride 10 --force-stride 5 --log-stride 15 --ramp-steps 10 --obstacle sdf --sdf-file $(GEOM_DIR)/formula1.sdfpoly --sdf-scale 65 --u-in 0.025 --tau 0.62 --output $(RESULTS_DIR)/smoke_formula1 --case-name smoke_formula1

airfoil: $(TARGET)
	./$(TARGET) --nx 768 --ny 320 --steps 1440 --frame-stride 4 --force-stride 4 --log-stride 180 --ramp-steps 240 --obstacle airfoil --chord 240 --camber 0.02 --camber-pos 0.40 --thickness 0.12 --airfoil-samples 700 --angle-deg 6 --u-in 0.032 --tau 0.62 --output $(RESULTS_DIR)/airfoil --case-name naca2412_aoa6_ultra

sports-car: $(TARGET)
	./$(TARGET) --nx 768 --ny 300 --steps 1440 --frame-stride 4 --force-stride 4 --log-stride 180 --ramp-steps 240 --obstacle sdf --sdf-file $(GEOM_DIR)/sports_car.sdfpoly --sdf-scale 230 --angle-deg 0 --u-in 0.030 --tau 0.62 --output $(RESULTS_DIR)/sports_car --case-name sports_car_sdf

formula1: $(TARGET)
	./$(TARGET) --nx 840 --ny 300 --steps 1800 --frame-stride 3 --force-stride 6 --log-stride 180 --ramp-steps 280 --obstacle sdf --sdf-file $(GEOM_DIR)/formula1.sdfpoly --sdf-scale 280 --angle-deg 0 --u-in 0.032 --tau 0.62 --output $(RESULTS_DIR)/formula1 --case-name formula1_sdf_presentation

gpu-benchmark: $(CUDA_TARGET)
	./$(CUDA_TARGET) --nx 512 --ny 192 --steps 300 --output $(RESULTS_DIR)/gpu_benchmark

gpu-benchmark-f1: $(CUDA_TARGET)
	test -f $(RESULTS_DIR)/formula1/frames/frame_000000.bin || (echo "Run make formula1 first"; exit 1)
	./$(CUDA_TARGET) --nx 840 --ny 300 --steps 600 --u-in 0.032 --tau 0.62 --solid-frame $(RESULTS_DIR)/formula1/frames/frame_000000.bin --solid-field-index 5 --output $(RESULTS_DIR)/formula1

clean:
	rm -rf $(BUILD_DIR)

-include $(DEP)
