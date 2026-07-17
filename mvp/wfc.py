"""Wave Function Collapse sampler for 4x4 box-pleated crease patterns.

Samples a locally-valid box-pleated crease pattern on a 4x4 grid (a 5x5
vertex lattice) by assigning a degree-4 vertex template to each of the 9
interior vertices. The vocabulary is the 96 oriented degree-4 templates from
``templates.expanded_templates()``. Local validity (Maekawa, Kawasaki,
big-little-big) is guaranteed by the templates; the sampler enforces
edge compatibility between neighbours and forbids crossing diagonals within a
cell.
"""

from __future__ import annotations

import copy
import random
from typing import Dict, List, Optional, Set, Tuple

from box_pleating import BoxPleatingPattern, CreaseType, Point

from templates import expanded_templates

GRID_SIZE = 4
INTERIOR: List[Tuple[int, int]] = [
    (i, j) for i in range(1, GRID_SIZE) for j in range(1, GRID_SIZE)
]

# direction_deg -> (di, dj); x increases right (i), y increases up (j)
DIR_OFFSET: Dict[int, Tuple[int, int]] = {
    0: (1, 0),
    45: (1, 1),
    90: (0, 1),
    135: (-1, 1),
    180: (-1, 0),
    225: (-1, -1),
    270: (0, -1),
    315: (1, -1),
}

Vocab = List[Dict]
Options = Dict[Tuple[int, int], Set[str]]
Assignment = Dict[Tuple[int, int], Optional[str]]


def _build_vocab_index() -> Tuple[Vocab, Dict[str, Dict]]:
    vocab = expanded_templates()
    by_id = {t["id"]: t for t in vocab}
    return vocab, by_id


VOCAB, BY_ID = _build_vocab_index()
ALL_IDS: Set[str] = set(BY_ID)

# Axis-only subset, retained for experiments. The sampler uses ALL_IDS plus
# crossing-diagonal pruning.
_AXIS = {0, 90, 180, 270}
AXIS_IDS: Set[str] = {
    tid for tid, t in BY_ID.items() if set(t["directions"]) <= _AXIS
}


def neighbour(vertex: Tuple[int, int], direction_deg: int) -> Tuple[int, int]:
    """Return the lattice vertex one step from ``vertex`` in ``direction_deg``."""
    di, dj = DIR_OFFSET[direction_deg]
    return vertex[0] + di, vertex[1] + dj


def opposite_direction(direction_deg: int) -> int:
    """Return ``direction_deg + 180``, mod 360."""
    return (direction_deg + 180) % 360


def edge_mv(template: Dict, direction_deg: int) -> Optional[int]:
    """Return the MV (+1/-1) of the template's crease in ``direction_deg``.

    Returns None if the template has no crease in that direction.
    ``template["edges"]`` is a frozenset of (direction, mv) pairs.
    """
    for direction, mv in template["edges"]:
        if direction == direction_deg:
            return mv
    return None


def compatible(tu: Dict, tv: Dict, direction_u_to_v: int) -> bool:
    """True iff templates ``tu`` (at u) and ``tv`` (at v) agree on the shared edge.

    They agree when both fire the shared edge with the same MV, or neither
    fires it. ``direction_u_to_v`` is the direction from u toward v; from v's
    side the back-direction is ``opposite_direction(direction_u_to_v)``.
    """
    return edge_mv(tu, direction_u_to_v) == edge_mv(
        tv, opposite_direction(direction_u_to_v)
    )


def _fires(tid: str, direction_deg: int) -> bool:
    return edge_mv(BY_ID[tid], direction_deg) is not None


def _corner_can_avoid(
    options: Options, corner: Tuple[int, int], direction_deg: int
) -> bool:
    """True if ``corner`` need not fire ``direction_deg``.

    Boundary corners (not in options) never fire from a template, so they
    avoid every direction.
    """
    if corner not in options:
        return True
    return any(not _fires(tid, direction_deg) for tid in options[corner])


def _corner_can_fire(
    options: Options, corner: Tuple[int, int], direction_deg: int
) -> bool:
    if corner not in options:
        return False
    return any(_fires(tid, direction_deg) for tid in options[corner])


def prune_crossing_diagonals(
    options: Options, worklist: List[Tuple[int, int]]
) -> bool:
    """Forbid both diagonals in the same unit cell.

    When / and \\ cross, add_crease creates a half-integer vertex whose four
    rays are two pairs of identical MV (each diagonal keeps one type). That
    always yields |M-V| in {0, 4}, never 2, which is Maekawa-impossible. So at
    most one diagonal per cell is allowed.

    For a cell with lower-left (x, y):
        /  present if SW fires 45 or NE fires 225
        \\ present if SE fires 135 or NW fires 315

    If / is forced, drop any option that would add \\, and vice versa; if both
    are forced, that is a contradiction.

    Mutates options and appends changed interior corners to worklist. Returns
    False on contradiction.
    """
    changed = True
    while changed:
        changed = False
        for x in range(GRID_SIZE):
            for y in range(GRID_SIZE):
                sw, se = (x, y), (x + 1, y)
                nw, ne = (x, y + 1), (x + 1, y + 1)

                must_slash = (not _corner_can_avoid(options, sw, 45)) or (
                    not _corner_can_avoid(options, ne, 225)
                )
                must_backslash = (not _corner_can_avoid(options, se, 135)) or (
                    not _corner_can_avoid(options, nw, 315)
                )
                if must_slash and must_backslash:
                    return False

                prunes: List[Tuple[Tuple[int, int], int]] = []
                if must_slash:
                    prunes.extend([(se, 135), (nw, 315)])
                if must_backslash:
                    prunes.extend([(sw, 45), (ne, 225)])

                for corner, direction in prunes:
                    if corner not in options:
                        continue
                    new_opts = {
                        tid for tid in options[corner] if not _fires(tid, direction)
                    }
                    if not new_opts:
                        return False
                    if new_opts != options[corner]:
                        options[corner] = new_opts
                        worklist.append(corner)
                        changed = True
    return True


def initial_options() -> Options:
    """Each interior vertex starts with every template id."""
    return {v: set(ALL_IDS) for v in INTERIOR}


def pick_min_entropy(options: Options, rng: random.Random):
    """Return an open vertex (len > 1) with the fewest remaining options.

    Ties are broken uniformly at random. Requires at least one open vertex.
    """
    min_len = float("inf")
    candidates = []

    for v, ids in options.items():
        length = len(ids)
        if length > 1:
            if length < min_len:
                min_len = length
                candidates = [v]
            elif length == min_len:
                candidates.append(v)

    return rng.choice(candidates)


def propagate(
    options: Options,
    vertex: Tuple[int, int],
    tid: str,
) -> bool:
    """Collapse ``vertex`` to ``tid`` and shrink neighbour options in place.

    Runs edge-compatibility AC-3, then crossing-diagonal pruning, repeating
    until both are quiet. Returns False on contradiction.
    """
    options[vertex] = {tid}
    worklist = [vertex]
    while worklist:
        u = worklist.pop()
        for d in DIR_OFFSET.keys():
            v = neighbour(u, d)
            if v not in options:
                continue
            new_v = {
                vid
                for vid in options[v]
                if any(compatible(BY_ID[uid], BY_ID[vid], d) for uid in options[u])
            }
            if not new_v:
                return False
            if new_v != options[v]:
                options[v] = new_v
                worklist.append(v)
        if not prune_crossing_diagonals(options, worklist):
            return False
    return True


def materialize(assignment: Assignment) -> BoxPleatingPattern:
    """Build a BoxPleatingPattern from a complete interior-vertex assignment.

    Each interior edge is shared by two vertices; to avoid inserting it twice,
    an interior-interior edge is added only from the lex-smaller endpoint,
    while edges toward a boundary vertex are always added (the boundary vertex
    is not in the assignment).

    Requires every INTERIOR vertex to map to a non-None template id.
    """
    pattern = BoxPleatingPattern(grid_size=GRID_SIZE)
    for (i, j), tid in assignment.items():
        for (direction, mv) in BY_ID[tid]["edges"]:
            end = neighbour((i, j), direction)
            if end not in assignment or (i, j) < end:
                pattern.add_crease(
                    Point(i, j),
                    Point(end[0], end[1]),
                    CreaseType.MOUNTAIN if mv == 1 else CreaseType.VALLEY,
                )
    return pattern


def wfc_sample(rng: Optional[random.Random] = None) -> Optional[BoxPleatingPattern]:
    """Run one WFC sample.

    Returns the pattern on success, or None if the search exhausts all options
    without a complete assignment. On a failed propagation the options are
    restored from a snapshot before the next alternative is tried.
    """
    rng = rng or random.Random()
    options = initial_options()
    stack: List[Tuple[Tuple[int, int], List[str], Options]] = []

    while True:
        if all(len(options[v]) == 1 for v in INTERIOR):
            assignment = {v: next(iter(options[v])) for v in INTERIOR}
            pattern = materialize(assignment)
            ok, report = pattern.is_valid_pattern()
            if not ok:
                raise RuntimeError(
                    "WFC produced an invalid pattern; likely a bug in "
                    f"propagate/materialize/crossing prune: {report}"
                )
            return pattern

        if any(len(options[v]) > 1 for v in INTERIOR):
            v = pick_min_entropy(options, rng)
            alts = list(options[v])
            rng.shuffle(alts)
            tid = alts.pop()
            before = copy.deepcopy(options)
            stack.append((v, alts, before))
            if propagate(options, v, tid):
                continue

        while True:
            if not stack:
                return None
            v, alts, before = stack.pop()
            options = before
            ok = False
            while alts:
                tid = alts.pop()
                snapshot = copy.deepcopy(options)
                stack.append((v, alts, snapshot))
                if propagate(options, v, tid):
                    ok = True
                    break
                stack.pop()
                options = snapshot
            if ok:
                break


def _check_geometry() -> None:
    assert neighbour((2, 2), 0) == (3, 2)
    assert neighbour((2, 2), 45) == (3, 3)
    assert neighbour((2, 2), 180) == (1, 2)
    assert opposite_direction(0) == 180
    assert opposite_direction(45) == 225
    assert opposite_direction(270) == 90

    c = next(
        t
        for t in VOCAB
        if t["origin"] == "C1" and set(t["directions"]) == {0, 90, 180, 270}
    )
    assert edge_mv(c, 0) in (1, -1)
    assert edge_mv(c, 45) is None

    ok_pair = None
    for tu in VOCAB:
        mu = edge_mv(tu, 0)
        if mu is None:
            continue
        for tv in VOCAB:
            mv = edge_mv(tv, 180)
            if mv == mu:
                ok_pair = (tu, tv, mu)
                break
        if ok_pair:
            break
    assert ok_pair is not None, "vocab should contain a matching east/west pair"
    tu, tv, _ = ok_pair
    assert compatible(tu, tv, 0) is True

    tv_bad = next(t for t in VOCAB if edge_mv(t, 180) is None)
    assert compatible(tu, tv_bad, 0) is False

    print(f"geometry OK  (vocab size={len(VOCAB)}, interior={len(INTERIOR)})")


def _check_min_entropy() -> None:
    rng = random.Random(0)
    opts: Options = {
        (1, 1): set(list(ALL_IDS)[:10]),
        (1, 2): set(list(ALL_IDS)[:3]),
        (1, 3): set(list(ALL_IDS)[:20]),
        (2, 1): set(list(ALL_IDS)[:1]),
        (2, 2): set(list(ALL_IDS)[:5]),
        (2, 3): set(list(ALL_IDS)[:8]),
        (3, 1): set(list(ALL_IDS)[:4]),
        (3, 2): set(list(ALL_IDS)[:6]),
        (3, 3): set(list(ALL_IDS)[:7]),
    }
    assert pick_min_entropy(opts, rng) == (1, 2)

    opts[(1, 3)] = set(list(ALL_IDS)[:3])
    seen = {pick_min_entropy(opts, random.Random(s)) for s in range(50)}
    assert seen <= {(1, 2), (1, 3)}, seen
    assert len(seen) == 2, "rng should eventually hit both tied vertices"
    print("min-entropy OK")


def _check_propagation() -> None:
    c = next(
        t
        for t in VOCAB
        if t["origin"] == "C1" and set(t["directions"]) == {0, 90, 180, 270}
    )
    mu = edge_mv(c, 0)
    assert mu is not None

    opts = initial_options()
    assert propagate(opts, (1, 1), c["id"]) is True
    assert opts[(1, 1)] == {c["id"]}

    for vid in opts[(2, 1)]:
        assert edge_mv(BY_ID[vid], 180) == mu, (
            f"east neighbour kept {vid} which does not match MV {mu} on west"
        )
    assert len(opts[(2, 1)]) < len(ALL_IDS)

    opts2 = initial_options()
    no_west = {tid for tid in ALL_IDS if edge_mv(BY_ID[tid], 180) is None}
    assert no_west, "expected some templates without a west edge"
    opts2[(2, 1)] = no_west
    assert propagate(opts2, (1, 1), c["id"]) is False

    print("propagation OK")


def _check_materialize() -> None:
    tid = next(
        t["id"]
        for t in VOCAB
        if t["origin"] == "C1" and set(t["directions"]) == {0, 90, 180, 270}
    )
    assignment: Assignment = {v: tid for v in INTERIOR}
    pattern = materialize(assignment)
    assert pattern.grid_size == GRID_SIZE
    assert len(pattern.creases) > 0

    segs = set()
    for c in pattern.creases:
        a = (c.start.x, c.start.y)
        b = (c.end.x, c.end.y)
        key = tuple(sorted([a, b]))
        assert key not in segs, f"duplicate crease {key}"
        segs.add(key)

    for v in INTERIOR:
        n = sum(
            1
            for c in pattern.creases
            if (c.start.x, c.start.y) == v or (c.end.x, c.end.y) == v
        )
        assert n == 4, f"{v} has degree {n}, expected 4"
    print("materialize OK")


def _check_crossing_prune() -> None:
    opts = initial_options()
    slash_tid = next(
        tid
        for tid, t in BY_ID.items()
        if t["origin"].startswith("A") and _fires(tid, 45) and not _fires(tid, 135)
    )
    assert propagate(opts, (1, 1), slash_tid) is True
    for tid in opts[(2, 1)]:
        assert not _fires(tid, 135), f"SE kept backslash template {tid}"
    for tid in opts[(1, 2)]:
        assert not _fires(tid, 315), f"NW kept backslash template {tid}"
    print(f"crossing-prune OK  (used {slash_tid})")


def _check_sampling() -> None:
    n_trials = 10
    successes = 0
    saw_diagonal = False
    for seed in range(n_trials):
        pattern = wfc_sample(random.Random(seed))
        if pattern is None:
            continue
        successes += 1
        ok, report = pattern.is_valid_pattern()
        assert ok, f"seed={seed} produced invalid pattern: {report}"
        for c in pattern.creases:
            dx = abs(c.end.x - c.start.x)
            dy = abs(c.end.y - c.start.y)
            if dx > 0 and dy > 0:
                saw_diagonal = True
    assert successes >= 1, f"expected some successes in {n_trials} trials, got {successes}"
    assert saw_diagonal, "expected at least one diagonal crease across successful samples"
    print(
        f"sampling OK  ({successes}/{n_trials} seeds valid; "
        f"vocab={len(ALL_IDS)}; saw_diagonal={saw_diagonal})"
    )


if __name__ == "__main__":
    _check_geometry()
    _check_min_entropy()
    _check_propagation()
    _check_materialize()
    _check_crossing_prune()
    _check_sampling()
