from math import isclose
import pytest
from modelcraft.residues import is_buffer, volume, weight


@pytest.mark.parametrize(
    "name,expected",
    [
        ("GOL", True),
        ("ALA", False),
    ],
)
def test_is_buffer(name: str, expected: bool):
    assert is_buffer(name) == expected


def test_volume():
    assert volume("HOH") == 18


@pytest.mark.parametrize(
    "name,expected",
    [
        ("ALA", 89.09),
        ("ASP", 132.09),
        ("ASN", 132.12),
        ("ASX", 130.12),
        ("UNK", 103.12),
    ],
)
def test_weight(name: str, expected: float):
    assert isclose(weight(name), expected, abs_tol=0.01)
