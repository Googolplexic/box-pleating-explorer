"""Degree-4 box-pleated vertex templates.

Enumerates the 10 canonical degree-4 box-pleating vertex templates,
modulo rotation/reflection of the vertex, and provides
`expanded_templates()` which generates the full WFC vocabulary by
applying the 16-element dihedral symmetry group D8 (8 rotations *
2 reflections) and de-duplicating.

Derivation:

- 8 cardinal directions, each a multiple of 45 deg.
- Sectors are 45 deg multiples summing to 360.
- Kawasaki: alpha_1 + alpha_3 = alpha_2 + alpha_4 = 180 deg.
  Up to cyclic shift this yields 3 sector-shape classes:
    A: (45,  45, 135, 135)
    B: (45,  90, 135,  90)   # two 90 deg sectors must sit OPPOSITE
    C: (90,  90,  90,  90)
  Note: (45, 90, 90, 135) fails Kawasaki (45+90=135 != 180) and is
  not a valid class.
- Maekawa: |M - V| = 2  =>  3M+1V or 1M+3V.
- Big-little-big: a sector strictly smaller than both neighbours must be
  bounded by opposite-MV creases. Only Class B has such a "little"
  sector (the 45 deg), so it constrains e1/e2 to opposite MV.
- Quotient by the symmetry of each shape:
    Class A has 2-fold rotational symmetry => 4 distinct MV templates.
    Class B has trivial symmetry           => 4 distinct MV templates.
    Class C has 4-fold rotational symmetry => 2 distinct MV templates.

Total: 4 + 4 + 2 = 10 templates.

Each template is a dict:
    {
        "id":         "A1",
        "shape":      "A",        # sector class
        "directions": [d0, d1, d2, d3],   # degrees, sorted increasing
        "mv":         [m0, m1, m2, m3],   # +1 mountain, -1 valley
        "notes":      "...",
    }

The directions are listed in counter-clockwise order around the vertex.
mv[i] is the MV label of the crease along directions[i].
"""

from __future__ import annotations

from typing import Dict, List

M = 1   # MOUNTAIN
V = -1  # VALLEY

TEMPLATES: List[Dict] = [
    # ---- Class A: sectors (45, 45, 135, 135) ----
    # Directions [0, 45, 90, 225] => sectors 45, 45, 135, 135 going CCW.
    {"id": "A1", "shape": "A", "directions": [0, 45, 90, 225],
     "mv": [V, M, M, M], "notes": "3M+1V, V on the e1=0deg crease"},
    {"id": "A2", "shape": "A", "directions": [0, 45, 90, 225],
     "mv": [M, V, M, M], "notes": "3M+1V, V on the e2=45deg crease"},
    {"id": "A3", "shape": "A", "directions": [0, 45, 90, 225],
     "mv": [M, V, V, V], "notes": "1M+3V, M on the e1=0deg crease"},
    {"id": "A4", "shape": "A", "directions": [0, 45, 90, 225],
     "mv": [V, M, V, V], "notes": "1M+3V, M on the e2=45deg crease"},

    # ---- Class B: sectors (45, 90, 135, 90) ----
    # Directions [0, 45, 135, 270] => sectors 45, 90, 135, 90 going CCW.
    # BLB: the 45deg sector is sector_1 (between e1=0 and e2=45);
    # both its neighbours are 90deg, so 45 is strictly smaller and
    # e1, e2 must have opposite MV. The 135deg sector is "big" (no BLB).
    {"id": "B1", "shape": "B", "directions": [0, 45, 135, 270],
     "mv": [V, M, M, M], "notes": "3M+1V, V on e1=0deg (e1!=e2 satisfied)"},
    {"id": "B2", "shape": "B", "directions": [0, 45, 135, 270],
     "mv": [M, V, M, M], "notes": "3M+1V, V on e2=45deg (e1!=e2 satisfied)"},
    {"id": "B3", "shape": "B", "directions": [0, 45, 135, 270],
     "mv": [M, V, V, V], "notes": "1M+3V, M on e1=0deg (e1!=e2 satisfied)"},
    {"id": "B4", "shape": "B", "directions": [0, 45, 135, 270],
     "mv": [V, M, V, V], "notes": "1M+3V, M on e2=45deg (e1!=e2 satisfied)"},

    # ---- Class C: sectors (90, 90, 90, 90) ----
    # Directions [0, 90, 180, 270] (orthogonal cross). The diagonal cross
    # [45, 135, 225, 315] is a 45-deg rotation of this and is treated as
    # equivalent here; the D8 expansion below recovers both orientations.
    {"id": "C1", "shape": "C", "directions": [0, 90, 180, 270],
     "mv": [V, M, M, M], "notes": "3M+1V orthogonal cross"},
    {"id": "C2", "shape": "C", "directions": [0, 90, 180, 270],
     "mv": [M, V, V, V], "notes": "1M+3V orthogonal cross"},
]


def template_count() -> Dict[str, int]:
    """Return per-class template counts (used for the assertion in
    the catalog renderer)."""
    counts: Dict[str, int] = {}
    for t in TEMPLATES:
        counts[t["shape"]] = counts.get(t["shape"], 0) + 1
    return counts


assert len(TEMPLATES) == 10, "Expected 10 templates (4+4+2)"
assert template_count() == {"A": 4, "B": 4, "C": 2}, template_count()


# ---------------------------------------------------------------------------
# D8 expansion: rotations + geometric reflection.
#
# The box-pleating direction set has dihedral symmetry of order 16 (8
# rotations by 45 deg, each composable with one reflection across the
# x-axis). Applying any of these 16 transforms to a valid template
# produces another valid template (Maekawa, Kawasaki, and BLB all
# preserve under rigid motions of the plane).
#
# Why no MV swap in this expansion: the 10 canonical templates already
# include both polarities (3M+1V and 1M+3V), so MV swap would only
# generate duplicates.
# ---------------------------------------------------------------------------

EdgeSet = frozenset  # frozenset[Tuple[int, int]]  -> set of (direction_deg, mv)


def _edges_of(template: Dict) -> EdgeSet:
    """Return a template's edges as a frozenset of (direction, mv) tuples."""
    return frozenset(zip(template["directions"], template["mv"]))


def _rotate(edges: EdgeSet, k_deg: int) -> EdgeSet:
    """Rotate every direction by k_deg degrees. MV labels unchanged."""
    return frozenset(((d + k_deg) % 360, mv) for d, mv in edges)


def _reflect(edges: EdgeSet) -> EdgeSet:
    """Reflect across the x-axis: d -> (-d) mod 360. MV labels unchanged."""
    return frozenset(((-d) % 360, mv) for d, mv in edges)


def expanded_templates() -> List[Dict]:
    """Generate the full WFC vocabulary: every D8 transform of every
    canonical template, de-duplicated.

    Each entry is a dict:
        {
            "id":         "A1_r0_f0",       # canonical_id _ rotation _ reflection
            "shape":      "A",
            "directions": [...sorted asc...],
            "mv":         [...same index as directions...],
            "edges":      frozenset of (direction, mv) tuples,
            "origin":     "A1",             # which canonical it came from
        }
    """
    seen: set = set()
    out: List[Dict] = []
    for canon in TEMPLATES:
        canon_edges = _edges_of(canon)
        for refl in (False, True):
            base = _reflect(canon_edges) if refl else canon_edges
            for k_idx in range(8):
                k_deg = k_idx * 45
                edges = _rotate(base, k_deg)
                if edges in seen:
                    continue
                seen.add(edges)
                sorted_edges = sorted(edges)
                out.append({
                    "id": f"{canon['id']}_r{k_idx}_f{int(refl)}",
                    "shape": canon["shape"],
                    "directions": [d for d, _ in sorted_edges],
                    "mv": [mv for _, mv in sorted_edges],
                    "edges": edges,
                    "origin": canon["id"],
                })
    return out


def expanded_count_by_shape() -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for t in expanded_templates():
        counts[t["shape"]] = counts.get(t["shape"], 0) + 1
    return counts


if __name__ == "__main__":
    exp = expanded_templates()
    print(f"Canonical templates: {len(TEMPLATES)}")
    print(f"Expanded templates : {len(exp)}")
    print(f"By shape           : {expanded_count_by_shape()}")
    print()
    print("First 6 expansions of A1:")
    for t in exp:
        if t["origin"] == "A1":
            print(f"  {t['id']}  dirs={t['directions']}  mv={t['mv']}")
