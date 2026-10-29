import gemmi
import numpy as np
import pytest

from modelcraft.mlgrow import next_residue, prev_residue
from modelcraft.tests.ccp4 import insulin_refmac

NAMES = ("N", "CA", "C")


def _residue(positions) -> gemmi.Residue:
    residue = gemmi.Residue()
    for name, position in zip(NAMES, positions):
        atom = gemmi.Atom()
        atom.name = name
        atom.pos = gemmi.Position(*position)
        residue.add_atom(atom)
    return residue


def _max_error(placed, target: gemmi.Residue) -> float:
    return max(
        np.linalg.norm(p - np.array(target[name][0].pos.tolist()))
        for p, name in zip(placed, NAMES)
    )


@pytest.mark.parametrize("grow", [next_residue, prev_residue])
def test_grown_residues_land_on_the_model(grow):
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    right = total = 0
    for chain in refmac.structure[0]:
        residues = [r for r in chain if r.find_atom("CA", "*") is not None]
        for r1, r2 in zip(residues, residues[1:]):
            source, target = (r1, r2) if grow is next_residue else (r2, r1)
            placed, error = grow(grid, source)
            assert 0 < error < 2
            right += _max_error(placed, target) < 1
            total += 1
    assert total > 40
    assert right / total > 0.75


@pytest.mark.parametrize("grow", [next_residue, prev_residue])
def test_grown_chains_follow_the_model(grow):
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    fractions = []
    for chain in refmac.structure[0]:
        residues = [r for r in chain if r.find_atom("CA", "*") is not None]
        if grow is prev_residue:
            residues.reverse()
        for i, start in enumerate(residues[:-1]):
            current, steps = start, 0
            for target in residues[i + 1 :]:
                placed, _ = grow(grid, current)
                if _max_error(placed, target) > 1.5:
                    break
                current = _residue(placed)
                steps += 1
            fractions.append(steps / (len(residues) - 1 - i))
    assert np.mean(fractions) > 0.35
