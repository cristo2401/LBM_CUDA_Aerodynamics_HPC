"""Check short CPU runs for complete, finite outputs (not physical validation)."""
import csv
import json
import math
from pathlib import Path
import struct
import sys


def check(directory: Path) -> None:
    metadata = json.loads((directory / "metadata.json").read_text())
    cells = metadata["nx"] * metadata["ny"]
    expected_bytes = cells * len(metadata["fields"]) * 4
    frames = sorted((directory / "frames").glob("frame_*.bin"))
    if not frames or len(frames) != metadata["frames_written"]:
        raise ValueError(f"{directory}: incomplete frame export")
    for frame in frames:
        raw = frame.read_bytes()
        if len(raw) != expected_bytes:
            raise ValueError(f"{frame}: wrong size")
        if not all(math.isfinite(value) for (value,) in struct.iter_unpack("=f", raw)):
            raise ValueError(f"{frame}: non-finite values")
    for name in ("diagnostics.csv", "forces.csv"):
        with (directory / name).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or int(rows[-1]["step"]) != metadata["steps"]:
            raise ValueError(f"{directory / name}: missing final step")
        if not all(math.isfinite(float(v)) for row in rows for v in row.values()):
            raise ValueError(f"{directory / name}: non-finite values")
        if name == "diagnostics.csv" and any(
            float(row["min_rho"]) <= 0 or int(row["finite"]) != 1 for row in rows
        ):
            raise ValueError(f"{directory}: invalid density/finite diagnostic")
    print(f"PASS: {directory} ({len(frames)} frames)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python3 scripts/check_smoke_outputs.py RESULTS_DIR [...]")
    for argument in sys.argv[1:]:
        check(Path(argument))
