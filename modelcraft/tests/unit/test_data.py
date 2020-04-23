from typing import List
import gemmi
from modelcraft.data import DataItem, cells_are_equal, _combine_data_items
from modelcraft.tests import data_path


def _test_data_item(mtz: gemmi.Mtz, labels: List[str]):
    item = DataItem(mtz, labels)
    assert cells_are_equal(item.cell, mtz.cell)
    assert item.spacegroup == mtz.spacegroup
    assert item.nreflections == mtz.nreflections
    assert len(item.columns) == len(labels) + 3
    assert item.label() == ",".join(labels)
    for i, label in enumerate(labels):
        assert item.label(i) == label


def test_1kb9_data_items():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    for labels in (["FREE"], ["FP", "SIGFP"], ["HLA", "HLB", "HLC", "HLD"]):
        _test_data_item(mtz, labels)


def test_hewl_data_items():
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    for labels in (
        ["F_New", "SIGF_New"],
        ["DANO_New", "SIGDANO_New"],
        ["F_New(+)", "SIGF_New(+)", "F_New(-)", "SIGF_New(-)"],
        ["IMEAN_New", "SIGIMEAN_New"],
        ["I_New(+)", "SIGI_New(+)", "I_New(-)", "SIGI_New(-)"],
        ["FreeR_flag"],
        ["FWT", "PHWT"],
        ["PHIB", "FOM"],
        ["HLA", "HLB", "HLC", "HLD"],
        ["HLanomA", "HLanomB", "HLanomC", "HLanomD"],
        ["parrot.ABCD.A", "parrot.ABCD.B", "parrot.ABCD.C", "parrot.ABCD.D"],
        ["parrot.F_phi.F", "parrot.F_phi.phi"],
    ):
        _test_data_item(mtz, labels)


def test_combine_data_items():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    types = [column.type for column in mtz.columns]
    fsigf = DataItem(mtz, ["FP", "SIGFP"])
    free = DataItem(mtz, ["FREE"])
    combined = _combine_data_items(fsigf, free)
    assert cells_are_equal(mtz.cell, combined.cell)
    assert mtz.spacegroup == combined.spacegroup
    assert mtz.nreflections == combined.nreflections
    assert combined.column_labels() == ["H", "K", "L", "FP", "SIGFP", "FREE"]
    types = [column.type for column in combined.columns]
    assert types == ["H", "H", "H", "F", "Q", "I"]


#     columns = [
#         "ISYM_New",
#         "FPFOM",
#     ]

# assert len(hklin.columns) == len(columns)
# assert all(col in hklin for col in columns)
# assert ",".join(columns) in hklin
# assert "LABEL_NOT_IN_FILE" not in hklin
# assert hklin.frees == {"FREE"}
# assert hklin.fsigfs == {"FP,SIGFP"}
# assert hklin.abcds == {"HLA,HLB,HLC,HLD"}
# assert hklin.phifoms == set()
# assert hklin.fphis == set()


# def test_hewl():
#     hklin = DataFile(data_path("hewl_data.mtz"))
#     assert hklin.spacegroup.hm == "P 43 21 2"
#     assert hklin.num_reflections == 10055
#     assert round(hklin.resolution_high, 2) == 1.86
#     assert round(hklin.resolution_low, 2) == 55.23

#     assert len(hklin.columns) == len(columns)
#     assert all(col in hklin for col in columns)
#     assert hklin.frees == {"FreeR_flag"}
#     assert hklin.fsigfs == {"F_New,SIGF_New"}
#     assert hklin.abcds == {
#         "HLA,HLB,HLC,HLD",
#         "HLanomA,HLanomB,HLanomC,HLanomD",
#         "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
#     }
#     assert hklin.phifoms == {
#         "PHWT,FPFOM",
#         "parrot.F_phi.phi,FPFOM",
#         "PHIB,FPFOM",
#         "PHIB,FOM",
#         "parrot.F_phi.phi,FOM",
#         "PHWT,FOM",
#     }
#     assert hklin.fphis == {
#         "FWT,PHWT",
#         "parrot.F_phi.F,parrot.F_phi.phi",
#     }
