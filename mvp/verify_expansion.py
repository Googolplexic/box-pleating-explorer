"""Verify every expanded template passes is_valid_pattern() and render a
mosaic of A1's D8 expansions.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt

from box_pleating import BoxPleatingPattern, CreaseType, Point
from box_pleating.viz import render_matplotlib

from templates import expanded_templates

GRID = 4
CENTER = (2, 2)


def _direction_to_boundary(direction_deg: int):
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
        if not pattern.add_crease(center, Point(end_x, end_y), crease_type):
            raise RuntimeError(
                f"Template {template['id']}: add_crease rejected the "
                f"{direction_deg}deg crease from {CENTER} to ({end_x},{end_y})."
            )
    return pattern


def main() -> None:
    expanded = expanded_templates()
    print(f"Verifying {len(expanded)} expanded templates...")

    by_shape = {"A": 0, "B": 0, "C": 0}
    failures = []
    for t in expanded:
        pattern = build_template_pattern(t)
        is_valid, report = pattern.is_valid_pattern()
        by_shape[t["shape"]] += 1
        if not is_valid:
            failures.append((t["id"], report))

    print(f"By shape: {by_shape}")
    print(f"Valid:    {len(expanded) - len(failures)} / {len(expanded)}")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):")
        for tid, report in failures[:5]:
            print(f"  {tid}: {report}")
        raise SystemExit(1)

    a1_expansions = [t for t in expanded if t["origin"] == "A1"]
    rows, cols = 3, 4
    fig, axes = plt.subplots(rows, cols, figsize=(12, 9))
    for ax, t in zip(axes.flat, a1_expansions):
        render_matplotlib(build_template_pattern(t), ax=ax, show_vertices=True)
        ax.set_title(t["id"], fontsize=9)
    for ax in axes.flat[len(a1_expansions):]:
        ax.axis("off")
    fig.suptitle(f"All {len(a1_expansions)} D8 expansions of canonical A1", fontsize=12)
    fig.tight_layout()
    out_path = Path(__file__).parent / "expansion_A1.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved A1 expansion mosaic to: {out_path}")


if __name__ == "__main__":
    main()
