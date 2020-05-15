import gemmi
from modelcraft.jobs import Refmac
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from modelcraft.tests import data_path


def test_1kv9():
    structure = read_structure(data_path("1kv9_model.pdb"))
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)
    assert refmac.fsigf.nreflections == mtz.nreflections
    assert refmac.abcd.nreflections == mtz.nreflections
    assert refmac.fphi_best.nreflections == mtz.nreflections
    assert refmac.fphi_diff.nreflections == mtz.nreflections
    assert refmac.fphi_calc.nreflections == mtz.nreflections
    assert refmac.rwork_change < 0
