import urllib
import gemmi
from modelcraft.jobs.refmac import RefmacXray, RefmacMapToMtz
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import ccp4_path, in_temp_directory


def test_1rxf():
    pdb_path = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    structure = read_structure(pdb_path)
    mtz_path = ccp4_path("examples", "data", "1rxf.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "F,SIGF")
    freer = DataItem(mtz, "FreeR_flag")
    refmac = RefmacXray(structure=structure, fsigf=fsigf, freer=freer, cycles=1).run()
    assert refmac.structure is not None
    assert refmac.abcd.nreflections == mtz.nreflections
    assert refmac.fphi_best.nreflections == mtz.nreflections
    assert refmac.fphi_diff.nreflections == mtz.nreflections
    assert refmac.fphi_calc.nreflections == mtz.nreflections
    assert 0 < refmac.rwork < 0.30
    assert 0 < refmac.rfree < 0.32
    assert 0.73 < refmac.fsc < 1
    assert refmac.data_completeness == 95.261
    assert refmac.resolution_high == 1.501


@in_temp_directory
def test_map_to_mtz():
    url = "https://ftp.wwpdb.org/pub/emdb/structures/EMD-3488/map/emd_3488.map.gz"
    urllib.request.urlretrieve(url, "emd_3488.map.gz")
    density = gemmi.read_ccp4_map("emd_3488.map.gz").grid
    mtztomap = RefmacMapToMtz(density, resolution=3.2).run()
    assert mtztomap.fphi is not None
