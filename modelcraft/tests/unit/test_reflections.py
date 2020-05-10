from typing import List
import gemmi
import pytest
from modelcraft.reflections import expand_label, DataItem, _combine_data_items
from modelcraft.tests import data_path


@pytest.mark.parametrize(
    "columns",
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
def test_valid_1kv9_free_columns(columns: str):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    free = DataItem(mtz, columns)
    assert len(free.columns) == 4


@pytest.mark.parametrize(
    "columns",
    [
        "",
        "*",
        "label",
        "dataset/label",
        "crystal/dataset/label",
        "project/crystal/dataset/label",
    ],
)
def test_invalid_1kv9_columns(columns: str):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    with pytest.raises(ValueError):
        DataItem(mtz, columns)


def test_1kv9_resolution():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    item = DataItem(mtz, "FP,SIGFP")
    assert item.resolution == pytest.approx(1.8, abs=0.001)


@pytest.mark.parametrize(
    "label,expanded",
    [
        ("HL", "HLA,HLB,HLC,HLD"),
        ("HLanom", "HLanomA,HLanomB,HLanomC,HLanomD"),
        ("HLABCD.F_sigF", "HLABCD.F_sigF.F,HLABCD.F_sigF.sigF"),
        ("parrot.ABCD", "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D"),
        ("parrot.F_phi", "parrot.F_phi.F,parrot.F_phi.phi"),
        ("FreeR_flag", "FreeR_flag"),
        ("prefix.HL", "prefix.HLA,prefix.HLB,prefix.HLC,prefix.HLD"),
        ("prefix.label", "prefix.label"),
        ("prefix_ABCD", "prefix_ABCD.A,prefix_ABCD.B,prefix_ABCD.C,prefix_ABCD.D"),
        ("prefix_HL", "prefix_HLA,prefix_HLB,prefix_HLC,prefix_HLD"),
        ("x.y.HL", "x.y.HLA,x.y.HLB,x.y.HLC,x.y.HLD"),
        ("prefix.FreeR_flag", "prefix.FreeR_flag"),
        ("prefix.I_sigI", "prefix.I_sigI.I,prefix.I_sigI.sigI"),
        ("prefix.phi_fom", "prefix.phi_fom.phi,prefix.phi_fom.fom"),
    ],
)
def test_expand_label(label: str, expanded: str):
    assert expand_label(label) == expanded


def test_1kv9_dataitem_init_types():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    item = DataItem(mtz, "FP,SIGFP")  # str
    assert item.label() == "FP,SIGFP"
    item = DataItem(mtz, mtz.columns[4:6])  # gemmi.MtzColumns
    assert item.label() == "FP,SIGFP"
    item = DataItem(mtz, [mtz.columns[4], mtz.columns[5]])  # List[gemmi.Mtz.Column]
    assert item.label() == "FP,SIGFP"


@pytest.mark.parametrize(
    "types,expected_labels",
    [("I", ["FREE"]), ("FQ", ["FP,SIGFP"]), ("AAAA", ["HLA,HLB,HLC,HLD"])],
)
def test_1kv9_item_search(types: str, expected_labels: List[str]):
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    labels = [item.label() for item in DataItem.search(mtz, types)]
    assert labels == expected_labels


@pytest.mark.parametrize(
    "columns",
    [
        "HLA,HLB,HLC,HLD",
        "HL",
        "New/HL",
        "New/HLA,HLB,HLC,HLD",
        "HLanomA,HLanomB,HLanomC,HLanomD",
        "HLanom",
        "New/HLanom",
        "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
        "parrot.ABCD",
        "New/parrot.ABCD",
    ],
)
def test_hewl_abcd_columns(columns: str):
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    abcd = DataItem(mtz, columns)
    assert abcd.types == "AAAA"
    assert len(abcd.columns) == 7


@pytest.mark.parametrize(
    "columns", ["New/HLA,HLB,HLC,HLD", "Old/HLA,HLB,HLC,HLD", "New/HL", "Old/HL"],
)
def test_valid_columns_for_mtz_with_multiple_datasets(columns):
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    mtz.add_dataset("Old")
    mtz.add_column("HLA", "A")
    mtz.add_column("HLB", "A")
    mtz.add_column("HLC", "A")
    mtz.add_column("HLD", "A")
    abcd = DataItem(mtz, columns)
    assert abcd.types == "AAAA"
    assert len(abcd.columns) == 7


@pytest.mark.parametrize(
    "columns", ["HLA,HLB,HLC,HLD", "HL"],
)
def test_invalid_columns_for_mtz_with_multiple_datasets(columns):
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    mtz.add_dataset("Old")
    mtz.add_column("HLA", "A")
    mtz.add_column("HLB", "A")
    mtz.add_column("HLC", "A")
    mtz.add_column("HLD", "A")
    with pytest.raises(ValueError):
        DataItem(mtz, columns)


@pytest.mark.parametrize(
    "types,expected_labels",
    [
        ("I", ["FreeR_flag"]),
        ("FQ", ["F_New,SIGF_New"]),
        ("FP", ["FWT,PHWT", "parrot.F_phi.F,parrot.F_phi.phi"]),
        ("PW", ["PHIB,FOM"]),
        (
            "AAAA",
            [
                "HLA,HLB,HLC,HLD",
                "HLanomA,HLanomB,HLanomC,HLanomD",
                "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
            ],
        ),
    ],
)
def test_hewl_item_search(types, expected_labels):
    mtz = gemmi.read_mtz_file(data_path("hewl_data.mtz"))
    labels = [item.label() for item in DataItem.search(mtz, types)]
    assert labels == expected_labels


def test_combine_data_items():
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    combined = _combine_data_items([fsigf, freer])
    assert combined.column_labels() == ["H", "K", "L", "FP", "SIGFP", "FREE"]
    assert mtz.nreflections == combined.nreflections
