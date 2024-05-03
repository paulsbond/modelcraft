from modelcraft.rszd import per_residue_rszd
from . import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    rszd = per_residue_rszd(refmac.structure, refmac.fphi_diff)
    assert len(rszd) > 0
    assert all(len(k) == 2 and v >= 0 for k, v in rszd.items())
