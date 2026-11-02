import gemmi
import numpy as np

from modelcraft.build import build, prepare_maps, to_structure
from modelcraft.build.chains import join, link, prune
from modelcraft.tests.ccp4 import insulin_contents, insulin_refmac, insulin_targets


def _chain_b():
    "N, CA and C of insulin chain B (30 residues)"
    chain = insulin_refmac().structure[0]["B"]
    return np.array(
        [
            [r[a][0].pos.tolist() for a in ("N", "CA", "C")]
            for r in chain
            if all(r.find_atom(a, "*") for a in ("N", "CA", "C"))
        ]
    )


def _grid():
    return prepare_maps(insulin_refmac().fphi_best)


def test_join_merges_overlapping_fragments():
    residues = _chain_b()
    cell = insulin_refmac().structure.cell
    spacegroup = gemmi.SpaceGroup("I 2 3")
    chains = join([residues[:18], residues[10:]], cell, spacegroup)
    assert len(chains) == 1
    assert np.abs(chains[0] - residues).max() < 1e-6


def test_prune_keeps_the_longer_of_two_overlapping_chains():
    residues = _chain_b()
    cell = insulin_refmac().structure.cell
    spacegroup = gemmi.SpaceGroup("I 2 3")
    chains = prune([residues, residues[5:15] + 0.5], cell, spacegroup)
    assert len(chains) == 1
    assert len(chains[0]) == len(residues)


def test_link_bridges_a_gap():
    grid, _ = _grid()
    residues = _chain_b()
    cell, spacegroup = grid.unit_cell, grid.spacegroup
    chains = link(grid, [residues[:12], residues[14:26]], cell, spacegroup)
    assert len(chains) == 1
    assert len(chains[0]) == 26
    assert np.linalg.norm(chains[0][:, 1] - residues[:26, 1], axis=1).max() < 1.5


def test_build_finds_most_of_insulin():
    grid, find_grid = _grid()
    main, _ = insulin_targets()
    sequences = [protein.sequence for protein in insulin_contents().proteins]
    chains = build(grid, find_grid, main, seed_pool=500, sequences=sequences)
    structure = insulin_refmac().structure
    structure.setup_cell_images()
    ns = gemmi.NeighborSearch(structure[0], structure.cell, 5).populate()
    built = np.concatenate(chains)
    correct = 0
    for residue in built:
        ca = gemmi.Position(*residue[1])
        marks = ns.find_atoms(ca, "\0", radius=1.0)
        correct += any(m.to_cra(structure[0]).atom.name == "CA" for m in marks)
    assert len(built) > 35
    assert correct / len(built) > 0.8


def test_structure_has_assigned_residue_names_and_numbers():
    residues = _chain_b()
    grid, _ = _grid()
    assignment = np.full((len(residues), 2), -1)
    assignment[2:5] = [[0, 0], [0, 1], [0, 3]]
    structure = to_structure(
        [residues], grid.unit_cell, grid.spacegroup, ["MKAW"], [assignment]
    )
    chain = structure[0][0]
    assert [r.name for r in chain][1:6] == ["UNK", "MET", "LYS", "TRP", "UNK"]
    assert [r.seqid.num for r in chain][:6] == [-1, 0, 1, 2, 4, 5]
