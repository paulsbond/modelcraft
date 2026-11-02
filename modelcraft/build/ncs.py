"""Copies of built chains in other copies of the molecule (NCS).

Chains assigned to the same sequence positions in different copies give an operator
that moves the chains near one copy onto the other, where they can fill in residues
that aren't built yet.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from ..geometry import superpose
from ..symmetry import Symmetry


def _pairs(chains, assignments, min_common, max_rmsd, min_shift):
    "(i, rotation, translation) superposing chain i onto another by their positions"
    positions = [
        {(int(s), int(p)): r for r, (s, p) in enumerate(a) if p >= 0}
        for a in assignments
    ]
    for i in range(len(chains)):
        for j in range(len(chains)):
            common = sorted(positions[i].keys() & positions[j].keys())
            if i == j or len(common) < min_common:
                continue
            a = chains[i][[positions[i][k] for k in common]].reshape(-1, 3)
            b = chains[j][[positions[j][k] for k in common]].reshape(-1, 3)
            if np.mean(np.linalg.norm(a - b, axis=1)) < min_shift:
                continue
            rotation, translation = superpose(a, b)
            moved = a @ rotation.T + translation
            rmsd = np.sqrt(np.mean(np.sum((moved - b) ** 2, axis=1)))
            if rmsd <= max_rmsd:
                yield i, rotation, translation


def copies(
    chains,
    assignments,
    cell,
    spacegroup,
    *,
    min_common=6,
    max_rmsd=1.0,
    min_shift=3.0,
    near=15.0,
):
    """New parts of the chains moved onto other copies of the molecule.

    A chain is moved onto another, with the chains within near of it, if they share
    at least min_common sequence positions, superpose within max_rmsd and are at
    least min_shift apart.
    """
    trees = [cKDTree(chain[:, 1]) for chain in chains]
    nearby = {}
    moved = []
    for i, rotation, translation in _pairs(
        chains, assignments, min_common, max_rmsd, min_shift
    ):
        if i not in nearby:
            nearby[i] = [
                k
                for k, chain in enumerate(chains)
                if np.isfinite(
                    trees[i].query(chain[:, 1], distance_upper_bound=near)[0]
                ).any()
            ]
        moved.extend(chains[k] @ rotation.T + translation for k in nearby[i])
    return new_parts(moved, chains, cell, spacegroup)


def new_parts(candidates, chains, cell, spacegroup, *, radius=2.0, overlap=3):
    """Runs of at least three candidate residues whose CAs aren't within radius of
    the chains or earlier parts, with overlap residues either side for joining."""
    if not candidates:
        return []
    sym = Symmetry(cell, spacegroup)
    if chains:
        sym.index(np.concatenate(chains)[:, 1], radius)
    added = np.zeros((0, 3))
    parts = []
    for candidate in candidates:
        cas = candidate[:, 1]
        built = sym.any_near(cas, radius) if chains else np.zeros(len(cas), bool)
        if len(added):
            built |= cdist(cas, added).min(axis=1) < radius
        start = None
        for r, new in enumerate(list(~built) + [False]):
            if new and start is None:
                start = r
            elif not new and start is not None:
                if r - start >= 3:
                    parts.append(candidate[max(start - overlap, 0) : r + overlap])
                    added = np.vstack([added, cas[start:r]])
                start = None
    return parts
