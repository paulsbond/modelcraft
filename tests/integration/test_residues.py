from modelcraft.residues import is_buffer, vdw_volume
from math import isclose


def test_is_buffer():
    assert is_buffer("GOL")
    assert not is_buffer("ALA")


def test_vdw_volume():
    assert isclose(vdw_volume("HOH"), 14.71, abs_tol=0.005)
    assert isclose(vdw_volume("HOH", include_hydrogrens=True), 29.19, abs_tol=0.005)
