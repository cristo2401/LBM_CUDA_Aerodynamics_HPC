"""Plot a saved solver field; see docs/PREVIEW.md for the reproducible run."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    directory = Path(sys.argv[1])
    destination = Path(sys.argv[2])
    metadata = json.loads((directory / "metadata.json").read_text())
    frame = sorted((directory / "frames").glob("frame_*.bin"))[-1]
    fields = np.fromfile(frame, dtype=np.float32).reshape(
        len(metadata["fields"]), metadata["ny"], metadata["nx"]
    )
    vorticity = fields[metadata["fields"].index("vorticity")]
    solid = fields[metadata["fields"].index("solid")] > 0.5
    visible = np.ma.masked_where(solid, vorticity)
    limit = max(float(np.max(np.abs(visible))), 1e-12)
    with plt.rc_context({"font.family": "DejaVu Sans", "svg.hashsalt": "lbm-preview"}):
        fig, ax = plt.subplots(figsize=(10, 4.8), layout="constrained")
        colors = plt.get_cmap("RdBu_r").copy()
        colors.set_bad("#111827")
        shown = ax.imshow(visible, origin="lower", cmap=colors,
                          vmin=-limit, vmax=limit, interpolation="nearest")
        ax.set(title=f"LBM D2Q9 · perfil alar 2D · paso {metadata['steps']}",
               xlabel="x (unidades de malla)", ylabel="y (unidades de malla)")
        fig.colorbar(shown, ax=ax, label="Vorticidad (unidades de malla)", shrink=0.85)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination, metadata={"Date": None})
        # A local PNG is convenient for reviewing the generated SVG.
        fig.savefig(destination.with_suffix(".png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python scripts/plot_preview.py RESULTS_DIR OUTPUT.svg")
    main()
