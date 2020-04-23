import gemmi
from modelcraft.data import DataItem
from modelcraft.sheetbend import Sheetbend
from modelcraft.tests import data_path
from modelcraft.model import model_stats


def test_1kv9():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, ["FP", "SIGFP"])
    free = DataItem(mtz, ["FREE"])
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    sheetbend = Sheetbend(fsigf, free, structure)
    assert model_stats(structure) == model_stats(sheetbend.structure)
