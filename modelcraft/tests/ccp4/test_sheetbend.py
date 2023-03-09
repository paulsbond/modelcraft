import gemmi
from modelcraft.jobs.sheetbend import Sheetbend
from modelcraft.reflections import DataItem
from modelcraft.structure import ModelStats, read_structure
from . import ccp4_path


def test_1rxf():
    mtz_path = ccp4_path("examples", "data", "1rxf.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "F,SIGF")
    freer = DataItem(mtz, "FreeR_flag")
    pdb_path = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    structure = read_structure(pdb_path)
    sheetbend = Sheetbend(structure=structure, fsigf=fsigf, freer=freer).run()
    assert ModelStats(structure) == ModelStats(sheetbend.structure)
