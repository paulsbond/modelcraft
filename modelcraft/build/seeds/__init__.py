"""Seed residues for building: search candidates scored by the grow network."""

import numpy as np

from ...stdres import STD_C, STD_CA, STD_N
from ..grow import next_residues, prev_residues


def _round_trips(grid, residues):
    "Round-trip predictions (back, forward) and a score (lower is better)"
    ahead, e1 = next_residues(grid, residues)
    back, e2 = prev_residues(grid, ahead)
    behind, e3 = prev_residues(grid, residues)
    forward, e4 = next_residues(grid, behind)
    distance = np.maximum(
        np.linalg.norm(back - residues, axis=2).max(axis=1),
        np.linalg.norm(forward - residues, axis=2).max(axis=1),
    )
    return back, forward, e1 + e2 + e3 + e4 + distance


def score_candidates(grid, candidates, pool=2500, refine=300):
    """Residues from the best search candidates, best first, scored by how confidently
    steps forward and back return to where they started."""
    order = np.argsort(candidates.zscores, kind="stable")[:pool]
    std = np.array([STD_N, STD_CA, STD_C])
    residues = np.matmul(std, candidates.rotations[order].transpose(0, 2, 1))
    residues += candidates.positions[order][:, None]
    _, _, score = _round_trips(grid, residues)
    residues = residues[np.argsort(score, kind="stable")[:refine]]
    for _ in range(3):
        back, forward, score = _round_trips(grid, residues)
        residues = (residues + back + forward) / 3
    _, _, score = _round_trips(grid, residues)
    return residues[np.argsort(score, kind="stable")]


def choose_seeds(residues, num_seeds, sym, model=None, clash=3.0):
    """The first residues with CAs (under symmetry) at least the clash radius from
    each other and from the model (a (n, 3, 3) array)."""
    has_model = model is not None and len(model) > 0
    if has_model:
        sym.index(model[:, 1], clash)
    seeds = []
    for residue in residues:
        if len(seeds) == num_seeds:
            break
        if has_model and sym.near(residue[1], clash):
            continue
        if any(sym.distance(residue[1], seed[1]) < clash for seed in seeds):
            continue
        seeds.append(residue)
    return seeds
