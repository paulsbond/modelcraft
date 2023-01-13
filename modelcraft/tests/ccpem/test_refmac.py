import gemmi
from modelcraft.jobs.refmac import RefmacMapToMtz
from . import density_path


def test_map_to_mtz():
    density = gemmi.read_ccp4_map(density_path())
    for blur in (0, 20, -20):
        mtztomap = RefmacMapToMtz(density, resolution=3.2, blur=blur).run()
        assert mtztomap.fphi is not None
