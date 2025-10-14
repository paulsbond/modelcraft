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
    contents.copies = None
    guessed = solvent_fraction(contents, cell, spacegroup, resolution=1.85)
    assert guessed == fraction


def test_1o5u():
    contents = AsuContents.from_pdbe("1o5u")
    assert "UNL" in {ligand.code for ligand in contents.ligands}
    cell = gemmi.UnitCell(39.627, 74.729, 42.21, 90, 90.02, 90)
    spacegroup = gemmi.SpaceGroup("P 1 21 1")
    fraction = solvent_fraction(contents, cell, spacegroup, resolution=1.83)
    assert fraction == pytest.approx(0.5, abs=0.1)
    contents.copies = None
    guessed = solvent_fraction(contents, cell, spacegroup, resolution=1.83)
    assert guessed == fraction


def test_9ckj():
    contents = AsuContents.from_pdbe("9ckj")
    assert "UNX" in set(contents.buffers)
    assert "A1AY8" in {ligand.code for ligand in contents.ligands}
    cell = gemmi.UnitCell(41.884, 76.025, 258.641, 90, 90, 90)
    spacegroup = gemmi.SpaceGroup("P 21 21 21")
    fraction = solvent_fraction(contents, cell, spacegroup, resolution=2.25)
    assert fraction == pytest.approx(0.5, abs=0.1)
    contents.copies = None
    guessed = solvent_fraction(contents, cell, spacegroup, resolution=2.25)
    assert guessed == fraction
