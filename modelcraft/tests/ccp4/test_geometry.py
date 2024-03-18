import gemmi
from modelcraft.structure import read_structure
from modelcraft.geometry import rmsz
from . import ccp4_path


def test_standard():
    path = ccp4_path("examples", "data", "2eck.pdb")
    structure = read_structure(path)
    assert 0 < rmsz(structure) < 2


def test_with_unl():
    path = ccp4_path("examples", "data", "2eck.pdb")
    structure = read_structure(path)
    residue = gemmi.Residue()
    residue.name = "UNL"
    residue.seqid = gemmi.SeqId(999, " ")
    structure[0][0].add_residue(residue)
    assert 0 < rmsz(structure) < 2
