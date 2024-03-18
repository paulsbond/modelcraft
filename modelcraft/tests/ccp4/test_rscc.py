from modelcraft.rscc import per_residue_rscc
from . import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    rscc = per_residue_rscc(refmac.structure, refmac.fphi_best)
    assert isinstance(rscc, dict)
    assert len(rscc) > 0
    assert all(isinstance(key, tuple) for key in rscc)
    assert all(isinstance(value, float) for value in rscc.values())
    assert all(-1 <= value <= 1 for value in rscc.values())
