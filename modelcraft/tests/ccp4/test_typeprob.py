import pytest

from modelcraft.tests.ccp4 import insulin_refmac
from modelcraft.typeprob import type_probabilities


def test_deposited_types_rank_highly():
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    top1 = top3 = total = 0
    for chain in refmac.structure[0]:
        for residue in chain:
            probabilities = type_probabilities(grid, residue)
            if residue.find_atom("CA", "*") is None:
                assert probabilities is None
                continue
            assert sum(probabilities.values()) == pytest.approx(1)
            ranked = sorted(probabilities, key=probabilities.get, reverse=True)
            top1 += ranked[0] == residue.name
            top3 += residue.name in ranked[:3]
            total += 1
    assert total > 40
    assert top1 / total > 0.35
    assert top3 / total > 0.65
