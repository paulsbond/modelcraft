import gemmi
import pytest
from modelcraft.reflections import find_column
from modelcraft.tests import data_path


# @pytest.mark.parametrize(
#     "label,split",
#     [
#         ("", (None, None, None, None)),
#         ("*", (None, None, None, None)),
#         ("label", (None, None, None, "label")),
#         ("*/label", (None, None, None, "label")),
#         ("/*/label", (None, None, None, "label")),
#         ("*/*/label", (None, None, None, "label")),
#         ("/*/*/label", (None, None, None, "label")),
#         ("*/*/*/label", (None, None, None, "label")),
#         ("/*/*/*/label", (None, None, None, "label")),
#         ("dataset/label", (None, None, "dataset", "label")),
#         ("/dataset/label", (None, None, "dataset", "label")),
#         ("crystal/dataset/label", (None, "crystal", "dataset", "label")),
#         ("/crystal/dataset/label", (None, "crystal", "dataset", "label")),
#         ("project/crystal/dataset/label", ("project", "crystal", "dataset", "label")),
#         ("/project/crystal/dataset/label", ("project", "crystal", "dataset", "label")),
#     ],
# )
# def test_split_column_label(label, split):
#     assert _split_column_label(label) == split


def test_find_column():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))

    def test_passing(pattern):
        assert find_column(mtz, pattern) == mtz.columns[3]

    def test_failing(pattern):
        with pytest.raises(ValueError):
            find_column(mtz, pattern)

    test_passing("FREE")
    test_passing("HKL_base/FREE")
    test_passing("*/FREE")
    test_passing("/*/FREE")
    test_passing("FREE")
    test_passing("*/FREE")
    test_passing("/*/FREE")
    test_passing("*/*/FREE")
    test_passing("/*/*/FREE")
    test_passing("*/*/*/FREE")
    test_passing("/*/*/*/*/*/*/*/*/*/*/*/*/FREE")
    test_passing("HKL_base/FREE")
    test_passing("HKL_base/HKL_base/FREE")
    test_passing("HKL_base/HKL_base/HKL_base/FREE")
    test_failing("")
    test_failing("*")
    test_failing("label")
    test_failing("dataset/label")
    test_failing("crystal/dataset/label")
    test_failing("project/crystal/dataset/label")


# def test_data_item_abc():
#     # with pytest.raises(TypeError):
#     DataItem()

# def _test_data_item(mtz: gemmi.Mtz, labels: List[str]):
#     item = DataItem(mtz, labels)
#     assert cells_are_equal(item.cell, mtz.cell)
#     assert item.spacegroup == mtz.spacegroup
#     assert item.nreflections == mtz.nreflections
#     assert len(item.columns) == len(labels) + 3
#     assert item.label() == ",".join(labels)
#     for i, label in enumerate(labels):
#         assert item.label(i) == label


# def test_1kb9_data_items():
#     mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
#     for labels in (["FREE"], ["FP", "SIGFP"], ["HLA", "HLB", "HLC", "HLD"]):
#         _test_data_item(mtz, labels)


# def test_hewl_data_items():
#     mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
#     for labels in (
#         ["F_New", "SIGF_New"],
#         ["DANO_New", "SIGDANO_New"],
#         ["F_New(+)", "SIGF_New(+)", "F_New(-)", "SIGF_New(-)"],
#         ["IMEAN_New", "SIGIMEAN_New"],
#         ["I_New(+)", "SIGI_New(+)", "I_New(-)", "SIGI_New(-)"],
#         ["FreeR_flag"],
#         ["FWT", "PHWT"],
#         ["PHIB", "FOM"],
#         ["HLA", "HLB", "HLC", "HLD"],
#         ["HLanomA", "HLanomB", "HLanomC", "HLanomD"],
#         ["parrot.ABCD.A", "parrot.ABCD.B", "parrot.ABCD.C", "parrot.ABCD.D"],
#         ["parrot.F_phi.F", "parrot.F_phi.phi"],
#     ):
#         _test_data_item(mtz, labels)


# def test_combine_data_items():
#     mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
#     types = [column.type for column in mtz.columns]
#     fsigf = DataItem(mtz, ["FP", "SIGFP"])
#     free = DataItem(mtz, ["FREE"])
#     combined = _combine_data_items(fsigf, free)
#     assert cells_are_equal(mtz.cell, combined.cell)
#     assert mtz.spacegroup == combined.spacegroup
#     assert mtz.nreflections == combined.nreflections
#     assert combined.column_labels() == ["H", "K", "L", "FP", "SIGFP", "FREE"]
#     types = [column.type for column in combined.columns]
#     assert types == ["H", "H", "H", "F", "Q", "I"]
