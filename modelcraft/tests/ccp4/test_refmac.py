import math
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
    refmac = Refmac(structure=structure, fsigf=fsigf, freer=freer, cycles=1).run()
    assert refmac.structure is not None
    assert math.isclose(refmac.initial_rwork, 0.3380, abs_tol=0.001)
    assert math.isclose(refmac.initial_rfree, 0.3429, abs_tol=0.001)
    assert math.isclose(refmac.initial_fsc, 0.9014, abs_tol=0.001)
    assert refmac.rwork < refmac.initial_rwork
    assert refmac.rfree < refmac.initial_rfree
    assert refmac.fsc > refmac.initial_fsc
    assert math.isclose(refmac.data_completeness, 95.261, abs_tol=0.01)
    assert math.isclose(refmac.resolution_high, 1.501, abs_tol=0.01)
