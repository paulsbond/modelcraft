from pytest import approx, mark
from modelcraft.monlib import atom_ids, in_library, is_buffer, volume, weight


@mark.parametrize(
    "code, expected",
    [
        ("HOH", {"O", "H1", "H2"}),
        ("GLY", {"N", "H", "H2", "H3", "CA", "HA3", "HA2", "C", "O", "OXT"}),
    ],
)
def test_atom_ids(code: str, expected: set):
    assert atom_ids(code) == expected


@mark.parametrize("code, expected", [("COM", True)])
def test_in_library(code: str, expected: bool):
    assert in_library(code) == expected


@mark.parametrize("code, expected", [("CL", True), ("ALA", False)])
def test_is_buffer(code: str, expected: bool):
    assert is_buffer(code) == expected


@mark.parametrize("code,expected", [("HOH", 18), ("2GP", 432)])
def test_volume(code: str, expected: float):
    assert volume(code) == approx(expected)


@mark.parametrize(
    "code,expected",
    [
        ("ALA", 89.09),
        ("ASP", 132.09),
        ("ASN", 132.12),
        ("ASX", 130.12),
        ("UNK", 103.12),
    ],
)
def test_weight(code: str, expected: float):
    assert weight(code) == approx(expected, abs=0.01)
