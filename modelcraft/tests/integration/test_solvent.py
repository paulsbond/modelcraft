import math
import pytest
from modelcraft.contents import Polymer, PolymerType
from modelcraft.solvent import _library_volume, _library_weight, _polymer_weight


@pytest.mark.parametrize(
    "code,expected",
    [
        ("ALA", 89.09),
        ("ASP", 132.09),
        ("ASN", 132.12),
        ("ASX", 130.12),
        ("UNK", 103.12),
    ],
)
def test_library_weight(code: str, expected: float):
    weight = _library_weight(code)
    assert math.isclose(weight, expected, abs_tol=0.01)


@pytest.mark.parametrize("code,expected", [("HOH", 18), ("2GP", 432)])
def test_library_volume(code: str, expected: float):
    volume = _library_volume(code)
    assert math.isclose(volume, expected, abs_tol=0.01)


def test_polymer_weight():
    polymer = Polymer("GG", polymer_type=PolymerType.PROTEIN)
    weight = _polymer_weight(polymer)
    assert math.isclose(weight, 132.12, abs_tol=0.01)
