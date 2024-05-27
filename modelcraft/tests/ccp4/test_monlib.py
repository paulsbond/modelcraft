import gemmi
from modelcraft.contents import PROTEIN_CODES, DNA_CODES, RNA_CODES
from modelcraft.monlib import atom_ids, in_library, group, is_protein, is_nucleic


def test_hoh_ids():
    ids = {"O", "H1", "H2"}
    assert atom_ids("HOH") == ids


def test_gly_ids():
    ids = {"N", "H", "H2", "H3", "CA", "HA3", "HA2", "C", "O", "OXT"}
    assert atom_ids("GLY") == ids


def test_com_ids():
    ids = {"C1", "C2", "S1", "S2", "O1S", "O2S", "O3S"}
    ids |= {"H11", "H12", "H21", "H22", "HS1", "HOS3"}
    assert atom_ids("COM") == ids


def test_in_library():
    for code in PROTEIN_CODES.values():
        assert in_library(code)
    for code in DNA_CODES.values():
        assert in_library(code)
    for code in RNA_CODES.values():
        assert in_library(code)
    assert not in_library("NOT_IN_MONLIB")


def test_group():
    assert group("GLY") == gemmi.ChemComp.Group.Peptide
    assert group("ALA") == gemmi.ChemComp.Group.Peptide
    assert group("MSE") == gemmi.ChemComp.Group.Peptide
    assert group("PRO") == gemmi.ChemComp.Group.PPeptide
    assert group("U") == gemmi.ChemComp.Group.Rna
    assert group("DT") == gemmi.ChemComp.Group.Dna
    assert group("HOH") == gemmi.ChemComp.Group.NonPolymer
    assert group("NOT_IN_MONLIB") is None


def test_protein():
    assert is_protein("GLY")
    assert is_protein("ALA")
    assert is_protein("MSE")
    assert is_protein("PRO")
    assert not is_protein("U")
    assert not is_protein("DT")
    assert not is_protein("HOH")
    assert not is_protein("NOT_IN_MONLIB")


def test_nucleic():
    assert not is_nucleic("GLY")
    assert not is_nucleic("ALA")
    assert not is_nucleic("PRO")
    assert is_nucleic("U")
    assert is_nucleic("DT")
    assert not is_nucleic("HOH")
    assert not is_nucleic("NOT_IN_MONLIB")
