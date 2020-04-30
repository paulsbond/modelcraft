import pytest
from modelcraft.contents import PolymerType, read_sequence_file
from modelcraft.tests import data_path


@pytest.mark.parametrize(
    "sequence,expected",
    [
        ("ABCDEFGHIKLMNPQRSTVWXYZ", PolymerType.PROTEIN),
        ("ACGU", PolymerType.RNA),
        ("ACGT", PolymerType.DNA),
        ("ACG", PolymerType.RNA),
        ("AAAA", PolymerType.PROTEIN),
        ("GGGG", PolymerType.PROTEIN),
    ],
)
def test_polymer_type_from_sequence(sequence: str, expected: PolymerType):
    assert PolymerType.from_sequence(sequence) == expected


def test_1kv9_read_sequence_file():
    path = data_path("1kv9_sequence.fasta")
    polymers = list(read_sequence_file(path))
    assert len(polymers) == 1
    assert polymers[0].label == "A"
    assert polymers[0].copies == 1
    assert polymers[0].polymer_type == PolymerType.PROTEIN
    assert len(polymers[0].sequence) == 668
    assert polymers[0].sequence[:5] == "AGVDE"
    assert polymers[0].sequence[-5:] == "HKAAP"


def test_hewl_read_sequence_file():
    path = data_path("hewl_sequence.fasta")
    polymers = list(read_sequence_file(path))
    assert len(polymers) == 1
    assert polymers[0].label == "hewl"
    assert polymers[0].copies == 1
    assert polymers[0].polymer_type == PolymerType.PROTEIN
    assert len(polymers[0].sequence) == 129
    assert polymers[0].sequence[:5] == "KVFGR"
    assert polymers[0].sequence[-5:] == "RGCRL"
