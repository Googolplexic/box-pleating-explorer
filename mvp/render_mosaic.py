"""Generate a 10x10 mosaic of WFC samples on a 4x4 grid.

Outputs mosaic.png next to this script.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from box_pleating.viz import render_matplotlib
from wfc import wfc_sample

N_SAMPLES = 100
COLS = 10
ROWS = N_SAMPLES // COLS
BASE_SEED = 0


def main() -> None:
    rng_seed = BASE_SEED
    patterns = []
    for i in range(N_SAMPLES):
        pattern = None
        for attempt in range(500):
            pattern = wfc_sample(__import__("random").Random(rng_seed + attempt))
            if pattern is not None:
                rng_seed += attempt + 1
                break
        if pattern is None:
            raise RuntimeError(f"failed to sample pattern {i} after 500 attempts")
        ok, report = pattern.is_valid_pattern()
        if not ok:
            raise RuntimeError(f"sample {i} invalid: {report}")
        patterns.append(pattern)

    fig, axes = plt.subplots(ROWS, COLS, figsize=(COLS * 1.4, ROWS * 1.4))
    for ax, idx, pattern in zip(axes.flat, range(N_SAMPLES), patterns):
        render_matplotlib(pattern, ax=ax)
        ax.set_title(f"#{idx}", fontsize=7)

    fig.suptitle(
        f"Locally valid box-pleated CPs — {N_SAMPLES}/{N_SAMPLES} samples "
        f"(4x4 grid, WFC)",
        fontsize=12,
    )
    fig.tight_layout()

    out_path = Path(__file__).parent / "mosaic.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
