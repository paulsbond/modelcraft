import gemmi
from modelcraft.jobs.refmac import Refmac
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import ccp4_path


def test_1rxf():
    pdb_path = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    structure = read_structure(pdb_path)
    mtz_path = ccp4_path("examples", "data", "1rxf.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "F,SIGF")
    freer = DataItem(mtz, "FreeR_flag")
    refmac = Refmac(structure, fsigf, freer, cycles=1).run()
    assert refmac.abcd.nreflections == mtz.nreflections
    assert refmac.fphi_best.nreflections == mtz.nreflections
    assert refmac.fphi_diff.nreflections == mtz.nreflections
    assert refmac.fphi_calc.nreflections == mtz.nreflections
    assert refmac.rwork < 100
