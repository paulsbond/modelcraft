import urllib
import gemmi
from modelcraft.jobs.emda import MapMask
from . import in_temp_directory


@in_temp_directory
def test_emda():
    url = "https://ftp.wwpdb.org/pub/emdb/structures/EMD-3488/map/emd_3488.map.gz"
    urllib.request.urlretrieve(url, "emd_3488.map.gz")
    density = gemmi.read_ccp4_map("emd_3488.map.gz").grid
    MapMask(density).run()
