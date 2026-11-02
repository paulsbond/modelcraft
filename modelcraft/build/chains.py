"""Growing, joining, pruning, linking and polishing chains of residues.

Residues are (3, 3) arrays of N, CA and C positions; chains are (n, 3, 3) arrays.
"""

import dataclasses

import numpy as np

from ..symmetry import Symmetry
from .grow import next_residues, prev_residues


@dataclasses.dataclass
class _Window:
    residues: np.ndarray  # (3, 3, 3): three residues of N, CA, C
    weight: float = 1.0
    alive: bool = True


def join(chains, cell, spacegroup, *, merge=2.0, link=2.0, min_windows=4):
    """Chains made by merging overlapping fragments through a graph of
    three-residue windows, as in Buccaneer's Ca_join."""
    sym = Symmetry(cell, spacegroup)
    windows = [
        _Window(chain[i : i + 3].copy())
        for chain in chains
        for i in range(len(chain) - 2)
    ]
    if not windows:
        return []
    centres = np.array([w.residues[1, 1] for w in windows])
    sym.index(centres, merge)
    for i, wi in enumerate(windows):
        if not wi.alive:
            continue
        for j, transform in sym.near(wi.residues[1, 1], merge):
            wj = windows[j]
            if j <= i or not wj.alive:
                continue
            moved = sym.apply(transform, wj.residues)
            if np.all(np.linalg.norm(moved[:, 1] - wi.residues[:, 1], axis=1) < merge):
                total = wi.weight + wj.weight
                wi.residues = (wi.weight * wi.residues + wj.weight * moved) / total
                wi.weight = total
                wj.alive = False
    alive = [w for w in windows if w.alive]
    centres = np.array([w.residues[1, 1] for w in alive])
    sym.index(centres, link)
    links = [[] for _ in alive]
    for i, wi in enumerate(alive):
        for j, transform in sym.near(wi.residues[2, 1], link):
            if j == i:
                continue
            moved = sym.apply(transform, alive[j].residues)
            if np.linalg.norm(moved[0, 1] - wi.residues[1, 1]) < link:
                links[i].append(j)
    weights = np.array([w.weight for w in alive])
    result = []
    used = np.zeros(len(alive), dtype=bool)
    while True:
        path = _best_path(links, used, weights)
        if len(path) < min_windows:
            break
        used[path] = True
        result.append(_assemble([alive[i].residues for i in path], sym))
    return result


def _best_path(links, used, weights):
    "Highest-scoring path through unused windows (score 1 + log weight per window)"
    n = len(links)
    if n == 0 or np.all(used):
        return []
    score = np.where(used, -np.inf, 1.0 + 0.1 * np.log(weights))
    best = score.copy()
    back = np.full(n, -1)
    order = _topological_order(links, used)
    for i in order:
        if used[i]:
            continue
        for j in links[i]:
            if not used[j] and best[i] + score[j] > best[j]:
                best[j] = best[i] + score[j]
                back[j] = i
    end = int(np.argmax(best))
    path = [end]
    seen = {end}
    while back[path[-1]] >= 0 and back[path[-1]] not in seen:
        path.append(back[path[-1]])
        seen.add(path[-1])
    return path[::-1]


def _topological_order(links, used):
    "Depth-first order of the windows, ignoring links that close cycles"
    n = len(links)
    state = np.zeros(n, dtype=int)
    order = []
    for start in range(n):
        if state[start] or used[start]:
            continue
        stack = [(start, iter(links[start]))]
        state[start] = 1
        while stack:
            node, children = stack[-1]
            for child in children:
                if state[child] == 0 and not used[child]:
                    state[child] = 1
                    stack.append((child, iter(links[child])))
                    break
            else:
                stack.pop()
                state[node] = 2
                order.append(node)
    return order[::-1]


def _assemble(windows, sym):
    "Chain of residues averaged over consecutive overlapping windows"
    placed = [windows[0]]
    for w in windows[1:]:
        placed.append(sym.nearest_copy(w, (0, 1), placed[-1][1, 1]))
    n = len(placed) + 2
    total = np.zeros((n, 3, 3))
    count = np.zeros(n)
    for i, w in enumerate(placed):
        total[i : i + 3] += w
        count[i : i + 3] += 1
    return total / count[:, None, None]


def prune(chains, cell, spacegroup, radius=2.0, min_length=6):
    """Chains cut where their CAs clash with those of longer chains, keeping pieces
    of at least min_length."""
    if not chains:
        return []
    sym = Symmetry(cell, spacegroup)
    owner = np.concatenate([np.full(len(c), i) for i, c in enumerate(chains)])
    index = np.concatenate([np.arange(len(c)) for c in chains])
    cas = np.concatenate([c[:, 1] for c in chains])
    sym.index(cas, radius)
    keep = [np.ones(len(c), dtype=bool) for c in chains]
    order = np.argsort([-len(c) for c in chains], kind="stable")
    rank = np.empty(len(chains), dtype=int)
    rank[order] = np.arange(len(chains))
    for p, (ci, ri) in enumerate(zip(owner, index)):
        if not keep[ci][ri]:
            continue
        for q, _ in sym.near(cas[p], radius):
            cj, rj = owner[q], index[q]
            if q == p or (cj == ci and abs(int(rj) - int(ri)) < 3):
                continue
            if keep[cj][rj] and rank[cj] > rank[ci]:
                keep[cj][rj] = False
            elif keep[cj][rj] and rank[cj] < rank[ci]:
                keep[ci][ri] = False
                break
    return [
        chain[start:end]
        for chain, mask in zip(chains, keep)
        for start, end in _runs(mask, min_length)
    ]


def _runs(mask, min_length):
    "(start, end) of each run of at least min_length True values"
    runs, start = [], None
    for i, value in enumerate(list(mask) + [False]):
        if value and start is None:
            start = i
        elif not value and start is not None:
            if i - start >= min_length:
                runs.append((start, i))
            start = None
    return runs


def _paths(grid, residues, grow, relaxed_error, max_steps):
    "Paths (b, steps, 3, 3) grown from residues, and how long each is valid (b,)"
    paths = np.zeros((len(residues), max_steps, 3, 3))
    lengths = np.zeros(len(residues), dtype=int)
    active = np.arange(len(residues))
    current = residues
    for step in range(max_steps):
        if len(active) == 0:
            break
        placed, errors = grow(grid, current)
        ok = errors <= relaxed_error
        active, current = active[ok], placed[ok]
        paths[active, step] = current
        lengths[active] = step + 1
    return paths, lengths


@dataclasses.dataclass
class _Link:
    length: int  # residues in the bridge
    source: int  # chain whose end is linked
    target: int  # chain whose start is linked
    bridge: np.ndarray  # (length, 3, 3) in the frame of the chain it grew from
    keep_to: int = None  # last residue kept from the source (backward links)
    keep_from: int = 0  # first residue kept from the target (forward links)
    forward: bool = True


def _find_links(grid, chains, sym, *, relaxed_error, max_steps, meet):
    "Candidate links, shortest bridge first"
    found = []
    for forward in (True, False):
        residues = np.array([c[-1] if forward else c[0] for c in chains])
        grow = next_residues if forward else prev_residues
        paths, lengths = _paths(grid, residues, grow, relaxed_error, max_steps)
        # Residues that a path can land on: the first or last three of each chain
        ends = [c[:3] if forward else c[-3:] for c in chains]
        targets = np.concatenate(ends)
        owner = np.concatenate([np.full(len(e), i) for i, e in enumerate(ends)])
        index = np.concatenate(
            [
                np.arange(len(e)) + (0 if forward else len(c) - len(e))
                for e, c in zip(ends, chains)
            ]
        )
        sym.index(targets[:, 1], meet)
        for a in range(len(chains)):
            for k in range(lengths[a]):
                hit = None
                for t, transform in sym.near(paths[a, k, 1], meet):
                    moved = sym.apply(transform, targets[t])
                    if (
                        owner[t] != a
                        and np.linalg.norm(moved - paths[a, k], axis=1).max() < meet
                    ):
                        hit = t
                        break
                if hit is None:
                    continue
                b = owner[hit]
                if forward:
                    found.append(_Link(k, a, b, paths[a, :k], keep_from=index[hit]))
                else:
                    bridge = paths[a, :k][::-1]
                    found.append(
                        _Link(k, b, a, bridge, keep_to=index[hit], forward=False)
                    )
                break
    found.sort(key=lambda link: link.length)
    return found


def link(grid, chains, cell, spacegroup, *, relaxed_error=2.0, max_steps=10, meet=0.7):
    """Chains linked where growing from one chain's end, with a relaxed error limit,
    lands within meet of another chain's first or last three residues."""
    if len(chains) < 2:
        return chains
    sym = Symmetry(cell, spacegroup)
    outgoing, incoming = {}, {}
    for candidate in _find_links(
        grid, chains, sym, relaxed_error=relaxed_error, max_steps=max_steps, meet=meet
    ):
        a, b = candidate.source, candidate.target
        if a in outgoing or b in incoming or _same_group(a, b, outgoing):
            continue
        outgoing[a] = candidate
        incoming[b] = candidate
    linked = []
    for first in range(len(chains)):
        if first in incoming:
            continue
        chain = chains[first]
        current = first
        while current in outgoing:
            candidate = outgoing[current]
            if candidate.keep_to is not None:
                cut = len(chains[current]) - 1 - candidate.keep_to
                chain = chain[: len(chain) - cut]
            target = chains[candidate.target][candidate.keep_from :]
            if candidate.forward:
                pieces = [candidate.bridge, target]
            else:
                pieces = [np.concatenate([candidate.bridge, target])]
            for piece in pieces:
                if len(piece):
                    moved = sym.nearest_copy(piece, (0, 1), chain[-1][1])
                    chain = np.concatenate([chain, moved])
            current = candidate.target
        linked.append(chain)
    return linked


def _same_group(a, b, outgoing):
    "Whether linking the end of chain a to the start of chain b would close a loop"
    node = b
    while node != a and node in outgoing:
        node = outgoing[node].target
    return node == a


def grow_fragments(grid, fragments, max_error=0.6, max_steps=50):
    """Fragments (lists of residues) grown at both ends until the predicted error is
    above max_error or the chain turns back on itself."""
    fragments = [list(f) for f in fragments]
    for forward in (True, False):
        end = -1 if forward else 0
        grow_end = next_residues if forward else prev_residues
        active = list(range(len(fragments)))
        for _ in range(max_steps):
            if not active:
                break
            ends = np.array([fragments[i][end] for i in active])
            placed, errors = grow_end(grid, ends)
            still = []
            for i, residue, error in zip(active, placed, errors):
                if error > max_error:
                    continue
                cas = np.array([r[1] for r in fragments[i]])
                if np.min(np.linalg.norm(cas - residue[1], axis=1)) < 2.0:
                    continue
                if forward:
                    fragments[i].append(residue)
                else:
                    fragments[i].insert(0, residue)
                still.append(i)
            active = still
    return [np.array(f) for f in fragments]


def _predictions(grid, chains):
    """The chains' residues concatenated (n, 3, 3), and the grow network's
    predictions (n, 3, 3) and errors (n,) of each from the residues before and after,
    with infinite errors where there is no neighbour"""
    residues = np.concatenate(chains) if chains else np.zeros((0, 3, 3))
    ends = np.cumsum([len(c) for c in chains])
    starts = ends - [len(c) for c in chains]
    predictions = []
    directions = ((next_residues, -1, starts), (prev_residues, 1, ends - 1))
    for grow, offset, edges in directions:
        has = np.ones(len(residues), bool)
        has[edges[edges < len(residues)]] = False
        targets = np.flatnonzero(has)
        placed, errors = grow(grid, residues[targets + offset])
        all_placed = residues.copy()
        all_errors = np.full(len(residues), np.inf)
        all_placed[targets] = placed
        all_errors[targets] = errors
        predictions.append((all_placed, all_errors))
    return residues, predictions, ends[:-1]


def polish(grid, chains, iterations=2):
    "Chains with each residue averaged with the grow network's predictions of it"
    for _ in range(iterations):
        if not chains:
            break
        residues, predictions, splits = _predictions(grid, chains)
        total = residues * 4.0
        weight = np.full(len(residues), 4.0)
        for placed, errors in predictions:
            w = np.where(np.isfinite(errors), 1 / np.maximum(errors, 0.1) ** 2, 0)
            total += placed * w[:, None, None]
            weight += w
        chains = np.split(total / weight[:, None, None], splits)
    return chains


def consistency(grid, chains):
    "Per-residue scores of agreement with the grow network (lower is better)"
    if not chains:
        return []
    residues, predictions, splits = _predictions(grid, chains)
    total = np.zeros(len(residues))
    count = np.zeros(len(residues))
    for placed, errors in predictions:
        has = np.isfinite(errors)
        distance = np.linalg.norm(placed[has] - residues[has], axis=2).max(axis=1)
        total[has] += errors[has] + distance
        count[has] += 1
    return np.split(total / np.maximum(count, 1), splits)


def trusted_fragments(grid, chains, max_score=0.8, min_length=3):
    """Runs of residues from a model (e.g. from molecular replacement) that the grow
    network agrees with, after polishing them towards its predictions."""
    chains = polish(grid, [c for c in chains if len(c) >= min_length], 2)
    return [
        chain[start:end]
        for chain, scores in zip(chains, consistency(grid, chains))
        for start, end in _runs(scores < max_score, min_length)
    ]
