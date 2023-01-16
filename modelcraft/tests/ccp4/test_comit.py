import gemmi
from modelcraft.jobs.comit import Comit
from modelcraft.reflections import DataItem
from . import ccp4_path


def test_comit_after_comit():
    mtz_path = ccp4_path("examples", "data", "gere.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "FPHASED,SIGFPHASED")
    fphi = DataItem(mtz, "FB,PHIB")
    comit1 = Comit(fsigf, fphi).run()
    comit2 = Comit(fsigf, comit1.fphi).run()
    assert comit2.abcd.nreflections == fsigf.nreflections
    assert comit2.fphi.nreflections == fsigf.nreflections
