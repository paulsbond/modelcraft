import gemmi

from modelcraft.jobs.servalcat import ServalcatFw
from modelcraft.reflections import DataItem

from . import ccp4_path


def ianom_test(path: str, label: str) -> None:
    mtz = gemmi.read_mtz_file(path)
    observations = DataItem(mtz, label)
    servalcat = ServalcatFw(observations=observations).run()
    assert servalcat.fmean.nreflections == observations.nreflections
    assert servalcat.fanom.nreflections == observations.nreflections
    assert servalcat.imean.nreflections == observations.nreflections


def imean_test(path: str, label: str) -> None:
    mtz = gemmi.read_mtz_file(path)
    observations = DataItem(mtz, label)
    servalcat = ServalcatFw(observations=observations).run()
    assert servalcat.fmean.nreflections == observations.nreflections
    assert servalcat.fanom is None
    assert servalcat.imean is None


def test_1rxf():
    path = ccp4_path("examples", "data", "1rxf.mtz")
    imean_test(path, "I,SIGI")


def test_ianom():
    path = ccp4_path("examples", "data", "1vr7_lr_i.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_deuterolysin():
    path = ccp4_path("examples", "data", "deuterolysin.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_gere():
    path = ccp4_path("examples", "data", "gere.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")
    ianom_test(path, "I(+),SIGI(+),I(-),SIGI(-)")


def test_insulin():
    path = ccp4_path("examples", "data", "insulin.mtz")
    imean_test(path, "IMEAN,SIGIMEAN")
