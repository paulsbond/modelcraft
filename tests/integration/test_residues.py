from modelcraft.residues import is_buffer


def test_is_buffer(code, expected):
    assert is_buffer("GOL")
    assert not is_buffer("ALA")
