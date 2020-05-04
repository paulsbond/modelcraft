import gemmi
from modelcraft.jobs import Refmac
from modelcraft.reflections import FsigF, FreeRFlag
from modelcraft.tests import data_path


def test_1kv9():
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)
    assert refmac.fsigf.nreflections == mtz.nreflections
    assert refmac.abcd.nreflections == mtz.nreflections
    assert refmac.fphi_best.nreflections == mtz.nreflections
    assert refmac.fphi_diff.nreflections == mtz.nreflections
    assert refmac.fphi_calc.nreflections == mtz.nreflections
    assert refmac.initial_rwork > refmac.final_rwork
