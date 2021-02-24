from modelcraft.residues import is_buffer, volume


def test_is_buffer():
    assert is_buffer("GOL")
    assert not is_buffer("ALA")


def test_volume():
    assert volume("HOH") == 18
