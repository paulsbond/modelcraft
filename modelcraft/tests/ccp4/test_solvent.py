import gemmi
import pytest
from modelcraft.contents import AsuContents
from modelcraft.solvent import solvent_fraction


def test_1o6a():
    contents = AsuContents.from_pdbe("1o6a")
    cell = gemmi.UnitCell(61.481, 61.481, 113.148, 90, 90, 90)
    spacegroup = gemmi.SpaceGroup("P 41 21 2")
    fraction = solvent_fraction(contents, cell, spacegroup, resolution=1.85)
    assert fraction == pytest.approx(0.5, abs=0.1)
    assert contents.copies == 2
    contents.copies = None
    assert solvent_fraction(contents, cell, spacegroup, resolution=1.85) == fraction
