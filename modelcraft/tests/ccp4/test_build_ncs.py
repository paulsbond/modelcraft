import gemmi
import numpy as np

from modelcraft.build.ncs import copies, new_parts
from modelcraft.tests.ccp4 import insulin_refmac


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


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


CELL = gemmi.UnitCell(200, 200, 200, 90, 90, 90)
P1 = gemmi.SpaceGroup("P 1")


def test_missing_residues_are_copied_from_another_copy():
    residues = _chain_b()
    rotation, translation = _rotation(1.0), np.array([80.0, 10.0, -5.0])
    other = residues @ rotation.T + translation
    first = np.column_stack([np.zeros(30, int), np.arange(30)])
    second = first[:15]
    moved = copies([residues, other[:15]], [first, second], CELL, P1)
    assert len(moved) == 1
    # The new residues 15-29 of the second copy, with three residues of overlap
    assert len(moved[0]) == 18
    assert np.abs(moved[0] - other[12:]).max() < 1e-6


def test_chains_with_different_positions_are_not_copied():
    residues = _chain_b()
    other = residues @ _rotation(1.0).T + 80
    first = np.column_stack([np.zeros(30, int), np.arange(30)])
    second = np.column_stack([np.zeros(15, int), np.arange(15) + 30])
    assert copies([residues, other[:15]], [first, second], CELL, P1) == []


def test_new_parts_skip_built_and_repeated_residues():
    residues = _chain_b()
    parts = new_parts([residues, residues], [residues[:20]], CELL, P1)
    assert len(parts) == 1
    assert np.abs(parts[0] - residues[17:]).max() < 1e-6
