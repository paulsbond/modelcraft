import gemmi
import pytest
from modelcraft.reflections import (
    _combine_data_items,
    find_column,
    DataItem,
    FreeRFlag,
    FsigF,
    FPhi,
    ABCD,
    PhiFom,
)
from modelcraft.tests import data_path


@pytest.mark.parametrize(
    "pattern",
    [
        "FREE",
        "*/FREE",
        "/*/FREE",
        "*/*/FREE",
        "/*/*/FREE",
        "*/*/*/FREE",
        "/*/*/*/*/*/*/*/*/*/*/*/*/FREE",
        "HKL_base/FREE",
        "HKL_base/HKL_base/FREE",
        "HKL_base/HKL_base/HKL_base/FREE",
    ],
)
def test_1kv9_free_valid_find_column(pattern):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    assert find_column(mtz, pattern) == mtz.columns[3]


@pytest.mark.parametrize(
    "pattern",
    [
        "",
        "*",
        "label",
        "dataset/label",
        "crystal/dataset/label",
        "project/crystal/dataset/label",
    ],
)
def test_1kv9_invalid_find_column(pattern):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    with pytest.raises(ValueError):
        find_column(mtz, pattern)


def test_dataitem_search_exception():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    with pytest.raises(TypeError):
        list(DataItem.search(mtz))


def test_1kv9_dataitem_init_types():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    item = DataItem(mtz, "FP,SIGFP")  # str
    assert item.label() == "FP,SIGFP"
    item = DataItem(mtz, mtz.columns[4:6])  # gemmi.MtzColumns
    assert item.label() == "FP,SIGFP"
    item = DataItem(mtz, [mtz.columns[4], mtz.columns[5]])  # List[gemmi.Mtz.Column]
    assert item.label() == "FP,SIGFP"


@pytest.mark.parametrize(
    "item_type,expected_labels",
    [(FreeRFlag, ["FREE"]), (FsigF, ["FP,SIGFP"]), (ABCD, ["HLA,HLB,HLC,HLD"])],
)
def test_1kv9_item_search(item_type, expected_labels):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    labels = [item.label() for item in item_type.search(mtz)]
    assert labels == expected_labels


@pytest.mark.parametrize(
    "item_type,expected_labels",
    [
        (FreeRFlag, ["FreeR_flag"]),
        (FsigF, ["F_New,SIGF_New"]),
        (FPhi, ["FWT,PHWT", "parrot.F_phi.F,parrot.F_phi.phi"]),
        (PhiFom, ["PHIB,FOM"]),
        (
            ABCD,
            [
                "HLA,HLB,HLC,HLD",
                "HLanomA,HLanomB,HLanomC,HLanomD",
                "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
            ],
        ),
    ],
)
def test_hewl_item_search(item_type, expected_labels):
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    labels = [item.label() for item in item_type.search(mtz)]
    assert labels == expected_labels


def test_combine_data_items():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    combined = _combine_data_items([fsigf, freer])
    assert combined.column_labels() == ["H", "K", "L", "FP", "SIGFP", "FREE"]
    assert mtz.nreflections == combined.nreflections
