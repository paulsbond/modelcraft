from modelcraft.jobs.servalcat import (
    ServalcatNemap,
    ServalcatTrim,
    ServalcatRefine,
    ServalcatFsc,
)
from modelcraft.maps import read_map
from modelcraft.structure import read_structure
from . import halfmap1_path, halfmap2_path, density_path, mask_path, structure_path


def test_servalcat_nemap():
    halfmap1 = read_map(halfmap1_path())
    halfmap2 = read_map(halfmap2_path())
    mask = read_map(mask_path())
    ServalcatNemap(halfmap1, halfmap2, 3.2, mask).run()


def test_servalcat_trim():
    halfmap1 = read_map(halfmap1_path())
    halfmap2 = read_map(halfmap2_path())
    density = read_map(density_path())
    mask = read_map(mask_path())
    trimmed = ServalcatTrim(mask, [density, halfmap1, halfmap2]).run()
    assert len(trimmed.maps) == 3


def test_servalcat_refine():
    structure = read_structure(structure_path())
    density = read_map(density_path())
    ServalcatRefine(structure, 3.2, density=density, cycles=1).run()


def test_servalcat_halfmap_refine():
    structure = read_structure(structure_path())
    halfmap1 = read_map(halfmap1_path())
    halfmap2 = read_map(halfmap2_path())
    ServalcatRefine(structure, 3.2, halfmap1, halfmap2, cycles=1).run()


def test_servalcat_fsc():
    structure = read_structure(structure_path())
    density = read_map(density_path())
    result = ServalcatFsc(structure, 3.2, density=density).run()
    assert result.fsc is not None
