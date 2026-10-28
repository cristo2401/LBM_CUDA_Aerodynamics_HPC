# Regenerar la figura del README

La figura representa un transitorio corto del solver CPU, no datos externos ni una
solución validada o convergida. Los parámetros están expresados en unidades de malla.

Desde la raíz del repositorio, compila con `make` y ejecuta:

```bash
OMP_NUM_THREADS=2 ./ejecuciones/build/lbm2d \
  --nx 192 --ny 96 --steps 600 --frame-stride 600 \
  --force-stride 100 --log-stride 600 --ramp-steps 100 \
  --obstacle airfoil --chord 60 --camber 0.02 --camber-pos 0.40 \
  --thickness 0.12 --airfoil-samples 200 --angle-deg 6 \
  --u-in 0.032 --tau 0.62 \
  --output resultados/readme_airfoil --case-name readme_airfoil

python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy matplotlib
python scripts/plot_preview.py resultados/readme_airfoil docs/assets/airfoil-vorticity.svg
```

El script lee el último frame y el esquema de `metadata.json`, enmascara el sólido
y usa una escala de color simétrica para la vorticidad. También genera un PNG local
para inspección. Solo el SVG se necesita para el README.
