import urllib
import gemmi
from modelcraft.jobs.emda import MapMask
from modelcraft.jobs.servalcat import ServalcatTrim
from modelcraft.structure import read_structure
from . import in_temp_directory, pdbe_download


@in_temp_directory
def test_emda():
    url = "https://ftp.wwpdb.org/pub/emdb/structures/EMD-3488/map/emd_3488.map.gz"
    urllib.request.urlretrieve(url, "emd_3488.map.gz")
    density = gemmi.read_ccp4_map("emd_3488.map.gz").grid
    mapmask = MapMask(density).run()
    pdbe_download("5ni1.cif")
    structure = read_structure("5ni1.cif")
    ServalcatTrim(density, mapmask.mask, structure).run()
