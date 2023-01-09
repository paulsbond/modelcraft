from os import environ
from pathlib import Path
import gemmi
from modelcraft.jobs.servalcat import ServalcatTrim, ServalcatNemap, ServalcatRefine
from modelcraft.structure import read_structure


def test_servalcat_trim():
    directory = Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data")
    density = gemmi.read_ccp4_map(str(directory / "map/mrc/emd_3488.map")).grid
    mask = gemmi.read_ccp4_map(str(directory / "map/mrc/emd_3488_mask.mrc")).grid
    structure = read_structure(str(directory / "pdb/5ni1.cif"))
    ServalcatTrim(density, mask, structure).run()


def test_servalcat_nemap():
    directory = Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data")
    halfmap1_path = str(directory / "map/mrc/3488_run_half1_class001_unfil.mrc")
    halfmap2_path = str(directory / "map/mrc/3488_run_half2_class001_unfil.mrc")
    halfmap1 = gemmi.read_ccp4_map(halfmap1_path)
    halfmap2 = gemmi.read_ccp4_map(halfmap2_path)
    mask = gemmi.read_ccp4_map(str(directory / "map/mrc/emd_3488_mask.mrc")).grid
    ServalcatNemap(halfmap1, halfmap2, mask, resolution=3.2).run()


def test_servalcat_refine():
    directory = Path(environ["CCPEM"], "lib/py2/ccpem/src/ccpem_core/test_data")
    structure = read_structure(str(directory / "pdb/5ni1.cif"))
    density = gemmi.read_ccp4_map(str(directory / "map/mrc/emd_3488.map")).grid
    refined = ServalcatRefine(structure, density, resolution=3.2, cycles=2).run()
    assert refined.fsc > 0.8
