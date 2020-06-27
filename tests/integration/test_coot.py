from modelcraft.jobs import Prune, FixSideChains
from modelcraft.structure import ModelStats
from tests.integration import remove_logs, insulin_refmac


def test_insulin_prune():
    refmac = insulin_refmac()
    refmac.structure.remove_alternative_conformations()
    prune = Prune(refmac.structure, refmac.fphi_best, refmac.fphi_diff)
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(prune.structure)
    assert stats_out.residues < stats_in.residues
    remove_logs()


def test_insulin_fix_side_chains():
    refmac = insulin_refmac()
    refmac.structure.remove_alternative_conformations()
    sidechains = FixSideChains(refmac.structure, refmac.fphi_best, refmac.fphi_diff)
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(sidechains.structure)
    assert stats_out.residues == stats_in.residues
    remove_logs()
