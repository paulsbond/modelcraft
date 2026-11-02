"""Build protein backbone into a map with the grow network.

Residues are (3, 3) arrays of N, CA and C positions; chains are (n, 3, 3) arrays.
"""

import gemmi
import numpy as np

from ..sequence import PROTEIN_CODES
from ..symmetry import Symmetry
from . import ncs
from .chains import grow_fragments, join, link, polish, prune, trusted_fragments
from .seeds import choose_seeds, find, score_candidates
from .sequence import assign, residue_numbers

CYCLES = 3


def to_structure(
    chains, cell, spacegroup, sequences=None, assignments=None
) -> gemmi.Structure:
    "A structure of the chains, named and numbered from any sequence assignments"
    structure = gemmi.Structure()
    structure.cell = cell
    structure.spacegroup_hm = spacegroup.hm
    model = gemmi.Model("1")
    for c, residues in enumerate(chains):
        chain = gemmi.Chain(_chain_name(c))
        assignment = np.full((len(residues), 2), -1)
        if assignments is not None:
            assignment = assignments[c]
        numbers = residue_numbers(assignment[:, 1])
        for i, ncac in enumerate(residues):
            residue = gemmi.Residue()
            residue.name = "UNK"
            if assignment[i, 1] >= 0:
                code = sequences[assignment[i, 0]][assignment[i, 1]].upper()
                residue.name = PROTEIN_CODES.get(code, "UNK")
            residue.seqid = gemmi.SeqId(*numbers[i])
            for name, position in zip(("N", "CA", "C"), ncac):
                atom = gemmi.Atom()
                atom.name = name
                atom.element = gemmi.Element(name[0])
                atom.pos = gemmi.Position(*position)
                atom.b_iso = 30.0
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    return structure


def _chain_name(index: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if index < len(letters):
        return letters[index]
    return letters[index // len(letters) % len(letters)] + letters[index % len(letters)]


def prepare_maps(fphi, find_resolution: float = 2.0, spacing: float = 0.5):
    """Maps for growing and, cut to find_resolution with the sampling seeds.find
    needs, for finding"""
    grid = fphi.map(spacing=spacing)
    data = np.array(fphi, copy=True)
    d = fphi.cell.calculate_d_array(data[:, :3].astype(int))
    mtz = gemmi.Mtz(with_base=True)
    mtz.cell, mtz.spacegroup = fphi.cell, fphi.spacegroup
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    mtz.set_data(data[d >= find_resolution])
    return grid, mtz.transform_f_phi_to_map("F", "PHI", sample_rate=3.0)


def build(
    grid: gemmi.FloatGrid,
    find_grid: gemmi.FloatGrid,
    target,
    *,
    known=None,
    sequences=None,
    seed_pool: int = 2500,
    log=None,
):
    """Backbone chains built into the maps from prepare_maps.

    Chains from a known model (e.g. from molecular replacement) are kept where the
    grow network agrees with them. If sequences are given, chains are copied to other
    copies of the molecule between cycles.
    """
    cell, spacegroup = grid.unit_cell, grid.spacegroup
    sym = Symmetry(cell, spacegroup)
    num_seeds = find.default_num_find(find_grid, per_100_residues=40)
    candidates = find.search(find_grid, target)
    chains = []
    if known:
        chains = trusted_fragments(grid, known)
        if log:
            log(f"{sum(len(c) for c in chains)} residues trusted from the known model")
    if chains:
        sym.index(np.concatenate(chains)[:, 1], 3.0)
        away = np.array([not sym.near(p, 3.0) for p in candidates.positions])
        candidates = find.Candidates(
            candidates.rotations[away],
            candidates.positions[away],
            candidates.zscores[away],
        )
    candidates = score_candidates(
        grid,
        candidates,
        pool=max(seed_pool, 20 * num_seeds),
        refine=max(300, 3 * num_seeds),
    )
    for cycle in range(CYCLES):
        model = np.concatenate(chains) if chains else None
        seeds = [[s] for s in choose_seeds(candidates, num_seeds, sym, model=model)]
        fragments = grow_fragments(grid, [list(c) for c in chains] + seeds)
        chains = join(fragments, cell, spacegroup)
        chains = prune(chains, cell, spacegroup)
        chains = link(grid, chains, cell, spacegroup)
        chains = polish(grid, chains, 1)
        copied = 0
        if sequences and cycle < CYCLES - 1:
            assignments = assign(grid, chains, sequences)
            copies = ncs.copies(chains, assignments, cell, spacegroup)
            copies = trusted_fragments(grid, copies, max_score=1.2)
            copied = sum(len(c) for c in copies)
            chains = chains + copies
        if log:
            lengths = [len(c) for c in chains]
            log(
                f"cycle {cycle + 1}: {len(seeds)} seeds, {sum(lengths)} residues "
                f"in {len(lengths)} chains, longest {max(lengths, default=0)}"
                + (f", {copied} copied by NCS" if copied else "")
            )
    return chains
