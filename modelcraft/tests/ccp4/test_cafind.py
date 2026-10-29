import gemmi
import numpy as np

from modelcraft.cafind import default_num_find, find
from modelcraft.llk import TYPES, frame, reference_targets
from modelcraft.tests.ccp4 import insulin_fsigf, insulin_refmac


def _targets_and_grid():
    refmac = insulin_refmac()
    main, sides = reference_targets(insulin_fsigf(), refmac.abcd)
    fphi = refmac.fphi_best
    grid = fphi.transform_f_phi_to_map(fphi.label(0), fphi.label(1), sample_rate=3)
    return refmac.structure, main, sides, grid


def test_main_chain_target_scores_the_model():
    structure, main, sides, grid = _targets_and_grid()
    assert set(sides) == set(TYPES)
    frames = [
        frame(*(residue[name][0].pos.tolist() for name in ("N", "CA", "C")))
        for chain in structure[0]
        for residue in chain
        if all(residue.find_atom(name, "*") for name in ("N", "CA", "C"))
    ]
    rotations, translations = zip(*frames)
    model = main.correl(grid, rotations, translations)
    shuffled = main.correl(grid, rotations[::-1], translations)
    assert len(model) > 40
    assert np.median(model) < -0.5
    assert np.median(shuffled) > np.median(model) + 0.3


def test_found_residues_are_on_the_model():
    structure, main, _, grid = _targets_and_grid()
    found = find(grid, main)
    assert len(found) == default_num_find(grid)
    assert [f.score for f in found] == sorted(f.score for f in found)
    structure.setup_cell_images()
    ns = gemmi.NeighborSearch(structure[0], structure.cell, 5).populate()
    distances = []
    for f in found:
        ca = gemmi.Position(*f.positions()[1])
        marks = ns.find_atoms(ca, "\0", radius=2.0)
        distances.append(
            min(
                (
                    structure.cell.find_nearest_pbc_image(ca, m.pos, 0).dist()
                    for m in marks
                    if m.to_cra(structure[0]).atom.name == "CA"
                ),
                default=np.inf,
            )
        )
    assert distances[0] < 1.0
    assert sum(d < 1.5 for d in distances) >= 3
