import gemmi
from modelcraft.jobs.ctruncate import CTruncate
from modelcraft.reflections import DataItem
from . import ccp4_path


def fanom_test(path: str, label: str) -> None:
    mtz = gemmi.read_mtz_file(path)
    observations = DataItem(mtz, label)
    ctruncate = CTruncate(observations=observations).run()
    assert ctruncate.fmean.nreflections == observations.nreflections
    assert ctruncate.fanom == observations
    assert ctruncate.imean is None
    assert ctruncate.ianom is None


def imean_test(path: str, label: str) -> None:
    mtz = gemmi.read_mtz_file(path)
    observations = DataItem(mtz, label)
    ctruncate = CTruncate(observations=observations).run()
    assert ctruncate.fmean.nreflections == observations.nreflections
    assert ctruncate.fanom is None
    assert ctruncate.imean == observations
    assert ctruncate.ianom is None


def ianom_test(path: str, label: str) -> None:
    mtz = gemmi.read_mtz_file(path)
    observations = DataItem(mtz, label)
    ctruncate = CTruncate(observations=observations).run()
    assert ctruncate.fmean.nreflections == observations.nreflections
    assert ctruncate.fanom.nreflections == observations.nreflections
    assert ctruncate.imean.nreflections == observations.nreflections
    assert ctruncate.ianom == observations


def test_1rxf():
    path = ccp4_path("examples", "data", "1rxf.mtz")
    imean_test(path, "I,SIGI")


def test_ianom():
    path = ccp4_path("examples", "data", "1vr7_lr_i.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_deuterolysin():
    path = ccp4_path("examples", "data", "deuterolysin.mtz")
    fanom_test(path, "F(+),SIGF(+),F(-),SIGF(-)")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_gere():
    path = ccp4_path("examples", "data", "gere.mtz")
    fanom_test(path, "F(+),SIGF(+),F(-),SIGF(-)")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_insulin():
    path = ccp4_path("examples", "data", "insulin.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")


def test_insulin_ssad():
    path = ccp4_path("examples", "data", "SSADinsulin.mtz")
    fanom_test(path, "F_CuKa(+),SIGF_CuKa(+),F_CuKa(-),SIGF_CuKa(-)")
