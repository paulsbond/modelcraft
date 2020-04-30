import gemmi
from modelcraft.reflections import FsigF, FreeRFlag
from modelcraft.sheetbend import Sheetbend
from modelcraft.structure import model_stats
from modelcraft.tests import data_path


def test_1kv9():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    sheetbend = Sheetbend(fsigf, freer, structure)
    assert model_stats(structure) == model_stats(sheetbend.structure)
