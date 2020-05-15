import gemmi
from modelcraft.reflections import DataItem
from modelcraft.jobs import Sheetbend
from modelcraft.structure import ModelStats, read_structure
from modelcraft.tests import data_path


def test_1kv9():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    structure = read_structure(data_path("1kv9_model.pdb"))
    sheetbend = Sheetbend(fsigf, freer, structure)
    assert ModelStats(structure) == ModelStats(sheetbend.structure)
