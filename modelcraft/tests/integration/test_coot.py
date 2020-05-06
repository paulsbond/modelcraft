import gemmi
from modelcraft.jobs import Refmac, Prune, FixSideChains
from modelcraft.tests import data_path
from modelcraft.structure import model_stats
from modelcraft.reflections import FsigF, FreeRFlag


def test_1kv9_prune():
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    stats_in = model_stats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    prune = Prune(structure, refmac.fphi_best, refmac.fphi_diff)
    stats_out = model_stats(prune.structure)
    assert stats_out.residues < stats_in.residues


def test_1kv9_fix_side_chains():
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    stats_in = model_stats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    sidechains = FixSideChains(structure, refmac.fphi_best, refmac.fphi_diff)
    stats_out = model_stats(sidechains.structure)
    assert stats_out.residues == stats_in.residues
