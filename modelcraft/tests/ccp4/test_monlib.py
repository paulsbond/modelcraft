import math

import gemmi
import pytest

from modelcraft.contents import DNA_CODES, PROTEIN_CODES, RNA_CODES
from modelcraft.monlib import MonLib


@pytest.fixture(name="monlib", scope="module")
def monlib_fixture():
    return MonLib({"COM", "2GP"}, include_standard=True)


def test_hoh_ids(monlib: MonLib):
    ids = {"O", "H1", "H2"}
    assert monlib.atom_ids("HOH") == ids


def test_gly_ids(monlib: MonLib):
    ids = {"N", "H", "H2", "H3", "CA", "HA3", "HA2", "C", "O", "OXT"}
    assert monlib.atom_ids("GLY") == ids


def test_com_ids(monlib: MonLib):
    ids = {"C1", "C2", "S1", "S2", "O1S", "O2S", "O3S"}
    ids |= {"H11", "H12", "H21", "H22", "HS1", "HOS3"}
    assert monlib.atom_ids("COM") == ids


def test_in_library(monlib: MonLib):
    for code in PROTEIN_CODES.values():
        assert code in monlib
    for code in DNA_CODES.values():
        assert code in monlib
    for code in RNA_CODES.values():
        assert code in monlib
    assert "NOT_IN_MONLIB" not in monlib


def test_group(monlib: MonLib):
    assert monlib.group("GLY") == gemmi.ChemComp.Group.Peptide
    assert monlib.group("ALA") == gemmi.ChemComp.Group.Peptide
    assert monlib.group("MSE") == gemmi.ChemComp.Group.Peptide
    assert monlib.group("PRO") == gemmi.ChemComp.Group.PPeptide
    assert monlib.group("U") == gemmi.ChemComp.Group.Rna
    assert monlib.group("DT") == gemmi.ChemComp.Group.Dna
    assert monlib.group("HOH") == gemmi.ChemComp.Group.NonPolymer
    assert monlib.group("NOT_IN_MONLIB") == gemmi.ChemComp.Group.Null


def test_protein(monlib: MonLib):
    assert monlib.is_protein("GLY")
    assert monlib.is_protein("ALA")
    assert monlib.is_protein("MSE")
    assert monlib.is_protein("PRO")
    assert not monlib.is_protein("U")
    assert not monlib.is_protein("DT")
    assert not monlib.is_protein("HOH")


def test_nucleic(monlib: MonLib):
    assert not monlib.is_nucleic("GLY")
    assert not monlib.is_nucleic("ALA")
    assert not monlib.is_nucleic("PRO")
    assert monlib.is_nucleic("U")
    assert monlib.is_nucleic("DT")
    assert not monlib.is_nucleic("HOH")


@pytest.mark.parametrize("code,expected", [("HOH", 18), ("2GP", 432)])
def test_volume(monlib: MonLib, code: str, expected: float):
    assert math.isclose(monlib.volume(code), expected, abs_tol=0.01)


@pytest.mark.parametrize(
    "code,expected",
    [
        ("ALA", 89.09),
        ("ASP", 132.09),
        ("ASN", 132.12),
        ("UNK", 103.12),
    ],
)
def test_weight(monlib: MonLib, code: str, expected: float):
    assert math.isclose(monlib.weight(code), expected, abs_tol=0.01)
