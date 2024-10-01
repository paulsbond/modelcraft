from ...structure import ModelStats
from ...prune import prune
from . import insulin_refmac


def test_insulin_prune():
    refmac = insulin_refmac()
    pruned = prune(
        structure=refmac.structure,
        fphi_best=refmac.fphi_best,
        fphi_diff=refmac.fphi_diff,
        fphi_calc=refmac.fphi_calc,
        threshold=-2,
    )
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(pruned)
    assert stats_out.residues < stats_in.residues
