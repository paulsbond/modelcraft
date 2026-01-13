from typing import List

import gemmi
import pytest

from modelcraft.reflections import (
    DataItem,
    column_refs,
    contract_label,
    convert_to_fsigf_and_phifom,
    expand_label,
)

from . import ccp4_path


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


@pytest.mark.parametrize(
    "label,contracted",
    [
        ("FreeR_flag", "FreeR_flag"),
        ("prefix.FreeR_flag", "prefix.FreeR_flag"),
        ("prefix.F_phi.F,prefix.F_phi.phi", "prefix.F_phi"),
        ("prefix.F_sigF.F,prefix.F_sigF.sigF", "prefix.F_sigF"),
        ("prefix.I_sigI.I,prefix.I_sigI.sigI", "prefix.I_sigI"),
        ("prefix.phi_fom.phi,prefix.phi_fom.fom", "prefix.phi_fom"),
        ("prefix.ABCD.A,prefix.ABCD.B,prefix.ABCD.C,prefix.ABCD.D", "prefix.ABCD"),
    ],
)
def test_contract_label(label: str, contracted: str):
    assert contract_label(label) == contracted


@pytest.mark.parametrize(
    "columns,project,crystal,dataset,label",
    [
        ("label", "", "", "", "label"),
        ("/label", "", "", "", "label"),
        ("/label/", "", "", "", "label"),
        ("[label]", "", "", "", "label"),
        ("[[label]]", "", "", "", "label"),
        ("/*/*/*/*/*/*/*/*/*/[[label]]//", "", "", "", "label"),
        ("dataset/label", "", "", "dataset", "label"),
        ("crystal/dataset/label", "", "crystal", "dataset", "label"),
        ("project/crystal/dataset/label", "project", "crystal", "dataset", "label"),
    ],
)
def test_single_column_ref(
    columns: str, project: str, crystal: str, dataset: str, label: str
):
    refs = column_refs(columns)
    assert refs[0].project == project
    assert refs[0].crystal == crystal
    assert refs[0].dataset == dataset
    assert refs[0].label == label


@pytest.mark.parametrize(
    "columns",
    [
        ("label1,label2"),
        ("[label1,label2]"),
        ("/*/*/*/label1,label2"),
        ("/*/*/*/[label1,label2]"),
    ],
)
def test_multiple_column_refs(columns: str):
    refs = column_refs(columns)
    for ref in refs:
        assert ref.project == ""
        assert ref.crystal == ""
        assert ref.dataset == ""
    assert refs[0].label == "label1"
    assert refs[1].label == "label2"


def test_column_refs_with_duplicate_labels():
    refs1 = column_refs("/xtal/peak/[F(+),SIGF(+),F(-),SIGF(-)]")
    refs2 = column_refs("/xtal/infl/[F(+),SIGF(+),F(-),SIGF(-)]")
    labels = ["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"]
    assert len(refs1) == len(labels)
    assert len(refs2) == len(labels)
    for i, label in enumerate(labels):
        assert refs1[i].project == ""
        assert refs2[i].project == ""
        assert refs1[i].crystal == "xtal"
        assert refs2[i].crystal == "xtal"
        assert refs1[i].dataset == "peak"
        assert refs2[i].dataset == "infl"
        assert refs1[i].label == label
        assert refs2[i].label == label


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


def test_convert_fphi_to_fsigf_and_phifom():
    path = ccp4_path("examples", "data", "gere.mtz")
    mtz = gemmi.read_mtz_file(path)
    fphi = DataItem(mtz, "FC,PHIC")
    fsigf, phifom = convert_to_fsigf_and_phifom(fphi)
    assert fsigf.label() == "F,SIGF"
    assert phifom.label() == "PHI,FOM"
