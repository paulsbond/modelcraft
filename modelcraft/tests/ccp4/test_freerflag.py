import gemmi
from modelcraft.jobs.freerflag import FreeRFlag
from . import ccp4_path


def test_1vr7():
    path = ccp4_path("examples", "data", "1vr7_lr_i.mtz")
    mtz = gemmi.read_mtz_file(path)
    result = FreeRFlag(mtz=mtz).run()
    assert result.freer is not None
