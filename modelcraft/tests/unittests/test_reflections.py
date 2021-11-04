import pytest
from modelcraft.reflections import column_refs, expand_label, contract_label


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
    for i in range(len(labels)):
        assert refs1[i].project == ""
        assert refs2[i].project == ""
        assert refs1[i].crystal == "xtal"
        assert refs2[i].crystal == "xtal"
        assert refs1[i].dataset == "peak"
        assert refs2[i].dataset == "infl"
        assert refs1[i].label == labels[i]
        assert refs2[i].label == labels[i]
