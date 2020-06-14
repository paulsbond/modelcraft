import os
import gemmi
from modelcraft.jobs import Refmac
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


def ccp4_path(*paths: str) -> str:
    if "CCP4" not in os.environ:
        raise EnvironmentError("CCP4 environment not set")
    return os.path.join(os.environ["CCP4"], *paths)


_insulin_refmac: Refmac = None


def insulin_refmac():
    global _insulin_refmac
    if _insulin_refmac is not None:
        return insulin_refmac
    mtz_path = ccp4_path("examples", "data", "insulin.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "F,SIGF")
    freer = DataItem(mtz, "FreeR_flag")
    pdb_path = ccp4_path("examples", "data", "insulin.pdb")
    structure = read_structure(pdb_path)
    refmac = Refmac(structure, fsigf, freer, cycles=0)
    _insulin_refmac = refmac
    return refmac
