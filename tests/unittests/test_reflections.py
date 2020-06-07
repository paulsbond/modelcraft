import pytest
from modelcraft.reflections import expand_label


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
