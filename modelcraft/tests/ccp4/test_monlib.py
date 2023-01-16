from modelcraft.monlib import atom_ids


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
