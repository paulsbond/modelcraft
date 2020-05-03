import gemmi
from modelcraft.jobs import Refmac, FindWaters
from modelcraft.tests import data_path
from modelcraft.structure import model_stats
from modelcraft.reflections import FsigF, FreeRFlag


def test_1kv9():
    structure = gemmi.read_structure(data_path("1kv9_model.pdb"))
    stats_in = model_stats(structure)
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    refmac = Refmac(structure, fsigf, freer, cycles=1)

    findwaters = FindWaters(structure, refmac.fphi_best)
    stats_out = model_stats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters > stats_in.waters
    assert stats_out.dummy_atoms == stats_in.dummy_atoms

    findwaters = FindWaters(structure, refmac.fphi_best, dummy=True)
    stats_out = model_stats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters == stats_in.waters
    assert stats_out.dummy_atoms > stats_in.dummy_atoms
