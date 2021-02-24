from modelcraft.polymer import Polymer, PolymerType
from math import isclose


def test_weight():
    polymer = Polymer("GG", polymer_type=PolymerType.PROTEIN)
    assert isclose(polymer.weight(), 132.12, abs_tol=0.01)
