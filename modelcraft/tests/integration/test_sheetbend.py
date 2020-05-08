import gemmi
from modelcraft.reflections import DataItem
from modelcraft.jobs import Sheetbend
from modelcraft.structure import model_stats
from modelcraft.tests import data_path


def test_1kv9():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    sheetbend = Sheetbend(fsigf, freer, structure)
    assert model_stats(structure) == model_stats(sheetbend.structure)
