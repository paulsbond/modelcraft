from modelcraft.contents import AsuContents, determine_polymer_type_from_sequence
from modelcraft.tests import data_path


def test_type_determination():
    assert determine_polymer_type_from_sequence("ABCDEFGHIKLMNPQRSTVWXYZ") == "protein"
    assert determine_polymer_type_from_sequence("ACGU") == "rna"
    assert determine_polymer_type_from_sequence("ACGT") == "dna"
    assert determine_polymer_type_from_sequence("ACG") == "rna"
    assert determine_polymer_type_from_sequence("AAAA") == "protein"
    assert determine_polymer_type_from_sequence("GGGG") == "protein"


def test_1kv9_sequence():
    contents = AsuContents()
    path = data_path("1kv9_sequence.fasta")
    contents.add_polymers_from_sequence_file(path)
    assert len(contents.polymers) == 1
    assert contents.polymers[0].polymer_type == "protein"
    assert len(contents.polymers[0].sequence) == 668
    assert contents.polymers[0].copies == "unknown"
