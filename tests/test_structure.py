from modelcraft.structure import _KNOWN_PROTEIN_RESIDUES


def test_number_of_known():
    assert len(_KNOWN_PROTEIN_RESIDUES) == 22


def test_unk_is_known():
    assert "UNK" in _KNOWN_PROTEIN_RESIDUES
