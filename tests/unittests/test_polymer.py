import pytest
from modelcraft.polymer import code1_to_code3, PolymerType


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
def test_polymer_type_guess_from_sequence(sequence: str, expected: PolymerType):
    assert PolymerType.guess_from_sequence(sequence) == expected


@pytest.mark.parametrize(
    "code1,polymer_type,expected",
    [
        ("A", PolymerType.PROTEIN, "ALA"),
        ("B", PolymerType.PROTEIN, "ASX"),
        ("J", PolymerType.PROTEIN, "UNK"),
        ("X", PolymerType.PROTEIN, "UNK"),
        ("Z", PolymerType.PROTEIN, "GLX"),
        ("A", PolymerType.RNA, "A"),
        ("J", PolymerType.RNA, "N"),
        ("X", PolymerType.RNA, "N"),
        ("A", PolymerType.DNA, "DA"),
        ("J", PolymerType.DNA, "DN"),
        ("X", PolymerType.DNA, "DN"),
    ],
)
def test_code1_to_code3(code1: str, polymer_type: PolymerType, expected: str):
    assert code1_to_code3(code1, polymer_type) == expected
