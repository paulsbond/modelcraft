from modelcraft.contents import PROTEIN_CODES, DNA_CODES, RNA_CODES
from modelcraft.monlib import atom_ids, in_library


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
