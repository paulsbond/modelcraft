import gemmi
import numpy as np
import pytest

from modelcraft.build import side
from modelcraft.build.sidechains import build, conformations
from modelcraft.tests.ccp4 import insulin_refmac

MAIN = ("N", "CA", "C")


def _main_chain(structure: gemmi.Structure) -> gemmi.Structure:
    "A copy of the protein residues with only N, CA and C"
    copy = gemmi.Structure()
    copy.cell = structure.cell
    copy.spacegroup_hm = structure.spacegroup_hm
    model = gemmi.Model("1")
    for chain in structure[0]:
        new_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if not all(residue.find_atom(a, "*") for a in MAIN):
                continue
            new_residue = gemmi.Residue()
            new_residue.name = residue.name
            new_residue.seqid = residue.seqid
            for name in MAIN:
                new_residue.add_atom(residue.find_atom(name, "*"))
            new_chain.add_residue(new_residue)
        model.add_chain(new_chain)
    copy.add_model(model)
    return copy


def test_lysine_conformations():
    n, ca, c = np.array([0.8431, 0, -1.204]), np.zeros(3), np.array([0.8775, 0, 1.253])
    atoms, options, frequencies = conformations("LYS", n, ca, c)
    assert atoms == ["CB", "CG", "CD", "CE", "NZ"]
    assert options.shape == (46, 5, 3)
    assert len(frequencies) == 46
    # CB moves little between rotamers
    cb = options[:, 0]
    assert np.all(np.linalg.norm(cb - cb[0], axis=-1) < 0.1)
    cb_ca = options[:, 0] - ca
    cb_cg = options[:, 0] - options[:, 1]
    cosines = np.sum(cb_ca * cb_cg, axis=1) / (
        np.linalg.norm(cb_ca, axis=1) * np.linalg.norm(cb_cg, axis=1)
    )
    assert np.all(np.abs(np.degrees(np.arccos(cosines)) - 114) < 30)


def test_insulin_side_chains_are_rebuilt():
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    structure = _main_chain(refmac.structure)
    structure[0][0][0].name = "UNK"
    build(structure, grid)
    deposited = {
        (chain.name, residue.seqid.num): residue for chain in refmac.structure[0]
        for residue in chain
    }  # fmt: skip
    close = total = 0
    for chain in structure[0]:
        for residue in chain:
            names = {a.name for a in residue}
            assert "O" in names
            ca = residue.find_atom("CA", "*")
            for atom in residue:
                if atom.name not in MAIN:
                    assert (atom.b_iso, atom.occ) == (ca.b_iso, ca.occ)
            if residue.name == "GLY":
                assert names == {"N", "CA", "C", "O"}
                continue
            assert "CB" in names
            reference = deposited[(chain.name, residue.seqid.num)]
            for atom in residue:
                expected = reference.find_atom(atom.name, "*")
                if expected and atom.name not in (*MAIN, "O", "CB"):
                    close += atom.pos.dist(expected.pos) < 1.0
                    total += 1
    assert structure[0][0][0].name == "UNK"
    assert {a.name for a in structure[0][0][0]} == {"N", "CA", "C", "O", "CB"}
    assert total > 100
    assert close / total > 0.6


def test_given_probabilities_match_and_side_chains_do_not_clash():
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    from_grid = _main_chain(refmac.structure)
    build(from_grid, grid)
    given = _main_chain(refmac.structure)
    ncac = np.array(
        [[r.find_atom(a, "*").pos.tolist() for a in MAIN]
         for chain in given[0] for r in chain]
    )  # fmt: skip
    build(given, probabilities=side.probabilities(grid, ncac))
    atoms = [
        (k, a.name, np.array(a.pos.tolist()))
        for k, r in enumerate(r for chain in given[0] for r in chain)
        for a in r
    ]
    expected = [a.pos.tolist() for chain in from_grid[0] for r in chain for a in r]
    assert np.allclose([p for *_, p in atoms], expected)
    positions = np.array([p for *_, p in atoms])
    owners = np.array([k for k, *_ in atoms])
    main = (*MAIN, "O", "CB")
    sides = [i for i, (_, name, _) in enumerate(atoms) if name not in main]
    d = np.linalg.norm(positions[sides][:, None] - positions[None], axis=-1)
    other = owners[sides][:, None] != owners[None]
    assert d[other].min() > 2.0


def test_top_rotamers_match_deposited_side_chains():
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    residues = [
        r for chain in refmac.structure[0] for r in chain
        if all(r.find_atom(a, "*") for a in MAIN)
    ]  # fmt: skip
    positions = np.array(
        [[r.find_atom(a, "*").pos.tolist() for a in MAIN] for r in residues]
    )
    joint = side.probabilities(grid, positions)
    close = total = 0
    for residue, row, ncac in zip(residues, joint, positions):
        atoms, options, _ = conformations(residue.name, *ncac)
        if not any(a != "CB" for a in atoms):
            continue
        probabilities = side.rotamer_probabilities(row[None], residue.name)[0]
        assert len(probabilities) == len(options)
        best = options[int(np.argmax(probabilities))]
        deposited = [residue.find_atom(a, "*") for a in atoms]
        if not all(deposited):
            continue
        # The atom after CB along the side chain (CG, OG, SG etc.)
        k = next(i for i, a in enumerate(atoms) if a != "CB")
        close += np.linalg.norm(best[k] - deposited[k].pos.tolist()) < 1.0
        total += 1
    assert total > 30
    assert close / total == pytest.approx(1, abs=0.15)
