import gemmi
from modelcraft.jobs.ctruncate import CTruncate
from modelcraft.reflections import DataItem
from tests.integration import ccp4_path


def test_imean():
    mtz_path = ccp4_path("examples", "data", "1rxf.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    observations = DataItem(mtz, "I,SIGI")
    ctruncate = CTruncate(observations=observations).run()
    assert ctruncate.fmean.nreflections == observations.nreflections
    assert ctruncate.fanom is None
    assert ctruncate.imean == observations
    assert ctruncate.ianom is None


def test_ianom():
    mtz_path = ccp4_path("examples", "data", "1vr7_lr_i.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    observations = DataItem(mtz, "I(+),SIGI(+),I(-),SIGI(-)")
    ctruncate = CTruncate(observations=observations).run()
    assert ctruncate.fmean.nreflections == observations.nreflections
    assert ctruncate.fanom.nreflections == observations.nreflections
    assert ctruncate.imean.nreflections == observations.nreflections
    assert ctruncate.ianom == observations
