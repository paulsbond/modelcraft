import gemmi
from modelcraft.jobs import Refmac, FindWaters
from modelcraft.tests import data_path
from modelcraft.structure import ModelStats, read_structure
from modelcraft.reflections import DataItem


def test_1kv9():
    structure = read_structure(data_path("1kv9_model.pdb"))
    stats_in = ModelStats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    findwaters = FindWaters(structure, refmac.fphi_best)
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters > stats_in.waters
    assert stats_out.dummy_atoms == stats_in.dummy_atoms

    findwaters = FindWaters(structure, refmac.fphi_best, dummy=True)
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters == stats_in.waters
    assert stats_out.dummy_atoms > stats_in.dummy_atoms
