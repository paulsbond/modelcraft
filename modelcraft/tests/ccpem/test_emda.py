from os import environ
from pathlib import Path
import gemmi
from modelcraft.jobs.emda import EmdaMapMask
from modelcraft.maps import read_map


def test_emda():
    directory = Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data")
    density = read_map(str(directory / "map/mrc/emd_3488.map"))
    mapmask = EmdaMapMask(density).run()
    assert isinstance(mapmask.mask, gemmi.Ccp4Map)
