# Proyecto LBM 2D

Simulador LBM D2Q9 en C++17 con casos de perfil alar, auto deportivo y Formula 1. La raiz queda ordenada en carpetas grandes:

- en librerias/cuda esta la implementación en CUDA, y en librerias/src/lbm_solver.cpp, esta la implementación en C++ CPU.


- `librerias/`: codigo C++, headers, kernels CUDA, benchmark CUDA y geometria importada.
- `resultados/`: frames binarios, `metadata.json`, fuerzas, diagnosticos y tiempos CUDA. Es local y se regenera con `make`.
- `ejecuciones/`: binarios y objetos de compilacion. Es local y se regenera con `make`.
- `notebooks/`: visualizador interactivo.
- `presentacion/`: previews, videos o material para exponer. Es local y no se sube al repo.
- `resto/`: espacio para material auxiliar que no sea codigo, resultados ni presentacion.

## Compilar

```bash
make
```

El binario principal queda en:

```text
ejecuciones/build/lbm2d
```

Sin OpenMP:

```bash
make OPENMP=0
```

Con CMake:

```bash
cmake -S . -B ejecuciones/cmake-build
cmake --build ejecuciones/cmake-build -j
```

## Ejecutar casos

Prueba corta:

```bash
make smoke
```

Perfil alar:

```bash
make airfoil
```

Vehiculos importados:

```bash
make vehicle-smoke
make sports-car
make formula1
```

Los contornos SDF/poligono estan en `librerias/geometries/`. Los resultados quedan separados por caso, por ejemplo:

- `resultados/airfoil`
- `resultados/formula1`
- `resultados/sports_car`
- `resultados/gpu_benchmark`

## Benchmark CUDA

El benchmark compara:

- `cpu_reference`
- `gpu_aos_naive`
- `gpu_soa_optimized`

Benchmark general:

```bash
make gpu-benchmark
```

Benchmark usando la mascara solida del Formula 1:

```bash
make formula1
make gpu-benchmark-f1
```

Exporta `gpu_kernel_timings.csv` y `gpu_kernel_comparison.json` dentro del directorio de resultados correspondiente.

## Salidas

Cada caso escribe:

- `metadata.json`: dimensiones, parametros, campos exportados y convencion de memoria.
- `frames/frame_XXXXXX.bin`: campos `float32` en orden `rho`, `ux`, `uy`, `pressure`, `vorticity`, `solid`, `phi`.
- `forces.csv`: `Fx`, `Fy`, `Cd`, `Cl`.
- `diagnostics.csv`: masa, densidades min/max, velocidad maxima y energia cinetica.

La indexacion es:

```text
index = x + nx * y
```

En Python cada campo se reconstruye con:

```python
field.reshape(ny, nx)
```

El notebook `notebooks/viewer.ipynb` lee por defecto `resultados/airfoil` y tiene una seccion especifica para `resultados/formula1` y la comparacion CUDA.
