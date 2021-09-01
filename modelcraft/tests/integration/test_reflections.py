from typing import List
import gemmi
from modelcraft.reflections import DataItem
from . import ccp4_path


def search_test(
    mtz: gemmi.Mtz, types: str, expected: List[str], sequential: bool = True
) -> None:
    items = DataItem.search(mtz, types, sequential)
    labels = [item.label() for item in items]
    assert labels == expected


def test_1rxf():
    path = ccp4_path("examples", "data", "1rxf.mtz")
    mtz = gemmi.read_mtz_file(path)
    search_test(mtz, "FQ", ["F,SIGF"])
    search_test(mtz, "JQ", ["I,SIGI"])
    search_test(mtz, "I", ["FreeR_flag"])


def test_1vr7():
    path = ccp4_path("examples", "data", "1vr7_lr_i.mtz")
    mtz = gemmi.read_mtz_file(path)
    search_test(mtz, "JQ", ["IMEAN,SIGIMEAN"])
    search_test(mtz, "KMKM", ["I(+),SIGI(+),I(-),SIGI(-)"])
    search_test(mtz, "I", [])


def test_deuterolysin():
    path = ccp4_path("examples", "data", "deuterolysin.mtz")
    mtz = gemmi.read_mtz_file(path)
    search_test(mtz, "FQ", ["F,SIGF"])
    search_test(mtz, "I", ["FreeR_flag"])
