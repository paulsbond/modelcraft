import functools
import os
import gemmi
from modelcraft.contents import AsuContents
from modelcraft.jobs import Refmac
from modelcraft.polymer import Polymer, PolymerType
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


def ccp4_path(*paths: str) -> str:
    return os.path.join(os.environ["CCP4"], *paths)


@functools.lru_cache(maxsize=None)
def insulin_fsigf():
    mtz_path = ccp4_path("examples", "data", "insulin.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    return DataItem(mtz, "F,SIGF")


@functools.lru_cache(maxsize=None)
def insulin_freer():
    mtz_path = ccp4_path("examples", "data", "insulin.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    return DataItem(mtz, "FreeR_flag")


@functools.lru_cache(maxsize=None)
def insulin_structure():
    pdb_path = ccp4_path("examples", "data", "insulin.pdb")
    return read_structure(pdb_path)


@functools.lru_cache(maxsize=None)
def insulin_refmac():
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    structure = insulin_structure()
    refmac = Refmac(structure, fsigf, freer, cycles=0)
    refmac.remove_files()
    return refmac


@functools.lru_cache(maxsize=None)
def insulin_contents():
    contents = AsuContents()
    a = Polymer(
        sequence="GIVEQCCASVCSLYQLENYCN",
        copies=1,
        polymer_type=PolymerType.PROTEIN,
    )
    b = Polymer(
        sequence="FVNQHLCGSHLVEALYLVCGERGFFYTPKA",
        copies=1,
        polymer_type=PolymerType.PROTEIN,
    )
    contents.polymers = [a, b]
    return contents
