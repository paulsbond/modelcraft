import pytest
from ...structure import read_structure
from . import ccp4_path


@pytest.mark.parametrize(
    "path",
    [
        ccp4_path("examples", "data", "insulin.pdb"),
        ccp4_path("examples", "toxd", "toxd.cif"),
    ],
)
def test_read_structure(path):
    structure = read_structure(path)
    assert len(structure) > 0
