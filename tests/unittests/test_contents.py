import pytest
from modelcraft.contents import (
    code1_to_code3,
    guess_sequence_type,
    PolymerType,
    sequences_in_file,
)


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
def test_guess_sequence_type(sequence: str, expected: PolymerType):
    assert guess_sequence_type(sequence) == expected


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


@pytest.mark.parametrize(
    "contents,expected",
    [
        ("AAAAAA", ["AAAAAA"]),
        ("AAAAAA\n", ["AAAAAA"]),
        (">\nAAAAAA", ["AAAAAA"]),
        (";\nAAAAAA", ["AAAAAA"]),
        (";ID\nAAAAAA", ["AAAAAA"]),
        ("\n>\nAAAAAA", ["AAAAAA"]),
        (">\nAAAAAA\n", ["AAAAAA"]),
        (">\nAAAAA A\n", ["AAAAAA"]),
        (">ID\nAAAAAA\n", ["AAAAAA"]),
        (">\n\nAAAAAA\n", ["AAAAAA"]),
        ("> ID\nAAAAAA\n", ["AAAAAA"]),
        (">\nAAAAAA---\n", ["AAAAAA"]),
        (">\n-AAA-AAA-\n", ["AAAAAA"]),
        (">\n\nAAAAAA\n\n", ["AAAAAA"]),
        (">P1;ID\n\nAAAAAA*", ["AAAAAA"]),
        (">\n\n\nAAAAAA\n\n", ["AAAAAA"]),
        (">P1;ID\n\nAAAAAA\n*", ["AAAAAA"]),
        (">\nAAAAAA\n;Comment\n", ["AAAAAA"]),
        (">P1;ID\nAAAAAA\nAAAAAA*", ["AAAAAA"]),
        (">P1;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">F1;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">D1;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">DL;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">DC;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">RL;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">RC;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">N3;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">N1;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">XX;ID\nComment\nAAAAAA*", ["AAAAAA"]),
        (">ID\nAAAAAA\nAAAAAA\n", ["AAAAAAAAAAAA"]),
        (">AA;ID\nAAAAAA\nAAAAAA*", ["AAAAAAAAAAAA"]),
        (">\nAAAAAA\n>\nAAAAAA\n", ["AAAAAA", "AAAAAA"]),
        (">P1;ID\nComment\nAAAAAA\n*\nComment", ["AAAAAA"]),
        (">P1;ID\nComment\nAAAAAA*\n\nComment", ["AAAAAA"]),
        (">ID\nAAAAAA\n>ID\nAAAAAA\n", ["AAAAAA", "AAAAAA"]),
        (";ID\n; Comment\nAAAAAA*\n\n>ID\nAAAAAA\n", ["AAAAAA", "AAAAAA"]),
        (">P1;ID\n\nAAAAAA\n*\n>P1;ID\n\nAAAAAA\n*\n", ["AAAAAA", "AAAAAA"]),
    ],
)
def test_sequences_in_file(contents: str, expected: list):
    assert sequences_in_file(contents) == expected
