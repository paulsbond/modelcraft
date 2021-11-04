from modelcraft.jobs.coot import Prune, FixSideChains, RsrMorph
from modelcraft.structure import ModelStats
from . import insulin_refmac


def test_insulin_prune():
    refmac = insulin_refmac()
    refmac.structure.remove_alternative_conformations()
    coot = Prune(
        structure=refmac.structure,
        fphi_best=refmac.fphi_best,
        fphi_diff=refmac.fphi_diff,
    ).run()
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(coot.structure)
    assert stats_out.residues < stats_in.residues


def test_insulin_fix_side_chains():
    refmac = insulin_refmac()
    refmac.structure.remove_alternative_conformations()
    coot = FixSideChains(
        structure=refmac.structure,
        fphi_best=refmac.fphi_best,
        fphi_diff=refmac.fphi_diff,
    ).run()
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(coot.structure)
    assert stats_out.residues == stats_in.residues


def test_insulin_rsr_morph():
    refmac = insulin_refmac()
    refmac.structure.remove_alternative_conformations()
    coot = RsrMorph(
        structure=refmac.structure,
        fphi_best=refmac.fphi_best,
    ).run()
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(coot.structure)
    assert stats_out.residues == stats_in.residues
