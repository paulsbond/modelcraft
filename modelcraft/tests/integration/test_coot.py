import gemmi
from modelcraft.jobs import Refmac, Prune, FixSideChains
from modelcraft.tests import data_path
from modelcraft.structure import ModelStats, read_structure
from modelcraft.reflections import DataItem


def test_1kv9_prune():
    structure = read_structure(data_path("1kv9_model.pdb"))
    stats_in = ModelStats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    prune = Prune(structure, refmac.fphi_best, refmac.fphi_diff)
    stats_out = ModelStats(prune.structure)
    assert stats_out.residues < stats_in.residues


def test_1kv9_fix_side_chains():
    structure = read_structure(data_path("1kv9_model.pdb"))
    stats_in = ModelStats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    sidechains = FixSideChains(structure, refmac.fphi_best, refmac.fphi_diff)
    stats_out = ModelStats(sidechains.structure)
    assert stats_out.residues == stats_in.residues
