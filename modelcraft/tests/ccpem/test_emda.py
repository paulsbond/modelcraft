from os import environ
from pathlib import Path
import gemmi
from modelcraft.jobs.emda import MapMask


def test_emda():
    directory = Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data")
    density = gemmi.read_ccp4_map(str(directory / "map/mrc/emd_3488.map"))
    mapmask = MapMask(density).run()
    assert isinstance(mapmask.mask, gemmi.Ccp4Mask)
