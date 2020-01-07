from modelcraft.contents import determine_polymer_type_from_sequence


def test_type_determination():
    assert determine_polymer_type_from_sequence("ABCDEFGHIKLMNPQRSTVWXYZ") == "protein"
    assert determine_polymer_type_from_sequence("ACGU") == "rna"
    assert determine_polymer_type_from_sequence("ACGT") == "dna"
    assert determine_polymer_type_from_sequence("ACG") == "rna"
    assert determine_polymer_type_from_sequence("AAAA") == "protein"
    assert determine_polymer_type_from_sequence("GGGG") == "protein"
