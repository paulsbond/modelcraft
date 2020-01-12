import modelcraft.residues as residues


def test_residue_count():
    assert len(residues.protein) == 20
    assert len(residues.rna) == 4
    assert len(residues.dna) == 4
