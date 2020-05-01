import pytest
from modelcraft.residues import PROTEIN, RNA, DNA


@pytest.mark.parametrize(
    "residues,expected", [(PROTEIN, 20), (RNA, 4), (DNA, 4)],
)
def test_residue_count(residues, expected):
    assert len(residues) == expected
