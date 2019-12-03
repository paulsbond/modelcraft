import modelcraft.gemmineer as gemmineer


def test_number_of_known():
    assert len(gemmineer._known_protein_residues) == 22


def test_unk_is_known():
    assert "UNK" in gemmineer._known_protein_residues
