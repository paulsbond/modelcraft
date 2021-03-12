import math
import pytest
from modelcraft.monomers import is_buffer, Monomers


_monomers = Monomers()


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
    assert _monomers.volume("HOH") == 18


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
    actual = _monomers.weight(name)
    assert math.isclose(actual, expected, abs_tol=0.01)
