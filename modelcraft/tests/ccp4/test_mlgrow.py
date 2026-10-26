import numpy as np
import pytest

from modelcraft.mlgrow import next_residue, prev_residue
from modelcraft.tests.ccp4 import insulin_refmac


@pytest.mark.parametrize("grow", [next_residue, prev_residue])
def test_grown_residues_land_on_the_model(grow):
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    right = total = 0
    for chain in refmac.structure[0]:
        residues = [r for r in chain if r.find_atom("CA", "*") is not None]
        for r1, r2 in zip(residues, residues[1:]):
            source, target = (r1, r2) if grow is next_residue else (r2, r1)
            placed, confidence = grow(grid, source)
            assert 0 < confidence < 2
            errors = [
                np.linalg.norm(p - np.array(target[name][0].pos.tolist()))
                for p, name in zip(placed, ("N", "CA", "C"))
            ]
            right += max(errors) < 1
            total += 1
    assert total > 40
    assert right / total > 0.80
