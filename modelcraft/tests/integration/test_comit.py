import gemmi
from modelcraft.reflections import DataItem
from modelcraft.jobs import Comit
from modelcraft.tests import data_path


def test_hewl():
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    fsigf = DataItem(mtz, "F_New,SIGF_New")
    fphi = DataItem(mtz, "FWT,PHWT")
    comit = Comit(fsigf, fphi)
    assert comit.abcd.nreflections == mtz.nreflections
    assert comit.fphi.nreflections == mtz.nreflections
