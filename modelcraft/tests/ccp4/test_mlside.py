import numpy as np

from modelcraft import mlside
from modelcraft.tests.ccp4 import insulin_refmac

MAIN = ("N", "CA", "C")


def _residues():
    "Insulin's deposited residues with N, CA and C, the map and their N, CA and C"
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    residues = [
        r for chain in refmac.structure[0] for r in chain
        if all(r.find_atom(a, "*") for a in MAIN)
    ]  # fmt: skip
    positions = np.array(
        [[r.find_atom(a, "*").pos.tolist() for a in MAIN] for r in residues]
    )
    return residues, grid, positions


def test_deposited_types_rank_highly():
    residues, grid, positions = _residues()
    joint = mlside.probabilities(grid, positions)
    assert np.allclose(joint.sum(axis=1), 1)
    types = mlside.type_probabilities(joint)
    assert types.shape == (len(residues), 20)
    assert np.allclose(types.sum(axis=1), 1)
    ranked = np.argsort(-types, axis=1)
    names = np.array(mlside.LABELS)[ranked]
    top1 = np.mean([n[0] == r.name for n, r in zip(names, residues)])
    top3 = np.mean([r.name in n[:3] for n, r in zip(names, residues)])
    assert len(residues) > 40
    assert top1 > 0.45
    assert top3 > 0.65


def test_rotamer_probabilities():
    residues, grid, positions = _residues()
    joint = mlside.probabilities(grid, positions)
    lys = mlside.rotamer_probabilities(joint, "LYS")
    assert lys.shape == (len(residues), 46)
    assert np.allclose(lys.sum(axis=1), 1)
    assert np.allclose(mlside.rotamer_probabilities(joint, "GLY"), 1)
    met = mlside.rotamer_probabilities(joint, "MET")
    assert np.allclose(mlside.rotamer_probabilities(joint, "MSE"), met)
    assert mlside.rotamer_probabilities(joint, "UNK") is None

