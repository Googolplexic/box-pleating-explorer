"""Verify every degree-4 template in templates.py and render a catalog mosaic.

For each template:
- Build a single-vertex BoxPleatingPattern centered at (2, 2) on a 4x4
  paper, with each crease running all the way to the paper boundary (so the
  only vertex requiring Maekawa+Kawasaki is the center).
- Assert is_valid_pattern() returns True.
- Render in a subplot.

Outputs catalog_degree4.png and prints a one-line report per template.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt

from box_pleating import BoxPleatingPattern, CreaseType, Point
from box_pleating.viz import render_matplotlib

from templates import TEMPLATES

GRID = 4
CENTER = (2, 2)


def _direction_to_boundary(direction_deg: int) -> Tuple[int, int]:
    """Project a direction outward from CENTER until it hits the paper
    boundary [0, GRID] x [0, GRID]. Box-pleating angles let us solve
    this exactly in integer coordinates on a 4x4 grid with center (2,2)."""
    cx, cy = CENTER
    dx = round(math.cos(math.radians(direction_deg)))
    dy = round(math.sin(math.radians(direction_deg)))
    steps_x = (GRID - cx) // dx if dx > 0 else (-cx // dx if dx < 0 else 10**9)
    steps_y = (GRID - cy) // dy if dy > 0 else (-cy // dy if dy < 0 else 10**9)
    steps = min(steps_x, steps_y)
    return cx + dx * steps, cy + dy * steps


def build_template_pattern(template) -> BoxPleatingPattern:
    pattern = BoxPleatingPattern(grid_size=GRID)
    center = Point(*CENTER)
    for direction_deg, mv in zip(template["directions"], template["mv"]):
        end_x, end_y = _direction_to_boundary(direction_deg)
        crease_type = CreaseType.MOUNTAIN if mv == 1 else CreaseType.VALLEY
        ok = pattern.add_crease(center, Point(end_x, end_y), crease_type)
        if not ok:
            raise RuntimeError(
                f"Template {template['id']}: add_crease rejected the "
                f"{direction_deg}deg crease from {CENTER} to ({end_x},{end_y})."
            )
    return pattern


def main() -> None:
    fig, axes = plt.subplots(2, 5, figsize=(15, 6))
    failures = []

    for ax, template in zip(axes.flat, TEMPLATES):
        pattern = build_template_pattern(template)
        is_valid, report = pattern.is_valid_pattern()

        render_matplotlib(pattern, ax=ax, show_vertices=True)
        title = f"{template['id']}: {template['shape']}  {'OK' if is_valid else 'FAIL'}"
        ax.set_title(title, fontsize=10)

        status = "OK" if is_valid else "FAIL"
        print(f"{template['id']:>3}  shape={template['shape']}  {status}  "
              f"directions={template['directions']}  mv={template['mv']}  "
              f"-- {template['notes']}")

        if not is_valid:
            failures.append((template["id"], report))

    fig.suptitle(
        "Degree-4 box-pleating vertex templates — "
        f"{len(TEMPLATES) - len(failures)}/{len(TEMPLATES)} valid",
        fontsize=12,
    )
    fig.tight_layout()
    out_path = Path(__file__).parent / "catalog_degree4.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved catalog to: {out_path}")

    if failures:
        print(f"\n{len(failures)} TEMPLATE(S) FAILED is_valid_pattern():")
        for tid, report in failures:
            print(f"  {tid}: {report}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
