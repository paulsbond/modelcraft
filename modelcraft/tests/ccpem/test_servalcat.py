import gemmi
from modelcraft.jobs.servalcat import ServalcatNemap, ServalcatTrim, ServalcatRefine
from modelcraft.structure import read_structure
from . import halfmap1_path, halfmap2_path, density_path, mask_path, structure_path


def test_servalcat_nemap():
    halfmap1 = gemmi.read_ccp4_map(halfmap1_path())
    halfmap2 = gemmi.read_ccp4_map(halfmap2_path())
    mask = gemmi.read_ccp4_mask(mask_path())
    ServalcatNemap(halfmap1, halfmap2, mask, 3.2).run()


def test_servalcat_trim():
    halfmap1 = gemmi.read_ccp4_map(halfmap1_path())
    halfmap2 = gemmi.read_ccp4_map(halfmap2_path())
    density = gemmi.read_ccp4_map(density_path())
    mask = gemmi.read_ccp4_mask(mask_path())
    structure = read_structure(structure_path())
    ServalcatTrim(mask, density, halfmap1, halfmap2, structure).run()


def test_servalcat_refine():
    structure = read_structure(structure_path())
    density = gemmi.read_ccp4_map(density_path())
    ServalcatRefine(structure, 3.2, density=density, cycles=1).run()


def test_servalcat_halfmap_refine():
    structure = read_structure(structure_path())
    halfmap1 = gemmi.read_ccp4_map(halfmap1_path())
    halfmap2 = gemmi.read_ccp4_map(halfmap2_path())
    ServalcatRefine(structure, 3.2, halfmap1, halfmap2, cycles=1).run()
