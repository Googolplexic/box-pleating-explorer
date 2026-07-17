"""Render the simplest valid degree-4 box-pleated vertex as a viz sanity check.

Builds a single interior vertex at (2, 2) on a 4x4 paper, with creases running
to the paper boundary at (4, 2), (2, 4), (0, 2), (2, 0). Boundary endpoints are
skipped by the verifier, so only the center must satisfy Maekawa + Kawasaki:

- MV pattern at the center: 3 Mountain + 1 Valley (Maekawa: |M - V| = 2).
- All sector angles are 90 deg (Kawasaki: 90+90 == 90+90 == 180).

Asserts the library considers it valid, then renders it to sanity_check.png.
"""

from pathlib import Path

import matplotlib.pyplot as plt

from box_pleating import BoxPleatingPattern, CreaseType, Point
from box_pleating.viz import render_matplotlib


def main() -> None:
    pattern = BoxPleatingPattern(grid_size=4)

    center = Point(2, 2)
    pattern.add_crease(center, Point(4, 2), CreaseType.MOUNTAIN)
    pattern.add_crease(center, Point(2, 4), CreaseType.VALLEY)
    pattern.add_crease(center, Point(0, 2), CreaseType.MOUNTAIN)
    pattern.add_crease(center, Point(2, 0), CreaseType.MOUNTAIN)

    is_valid, report = pattern.is_valid_pattern()
    assert is_valid, f"Expected valid pattern. Report: {report}"
    print("is_valid_pattern: True")

    out_path = Path(__file__).parent / "sanity_check.png"
    fig, ax = plt.subplots(figsize=(4, 4))
    render_matplotlib(pattern, ax=ax, show_vertices=True)
    ax.set_title("Degree-4 vertex (3M + 1V, 90° sectors)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved render to: {out_path}")


if __name__ == "__main__":
    main()
