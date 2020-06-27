from modelcraft.jobs import FindWaters
from modelcraft.structure import ModelStats
from tests.integration import remove_logs, insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    stats_in = ModelStats(refmac.structure)

    findwaters = FindWaters(refmac.structure, refmac.fphi_best)
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters > stats_in.waters
    assert stats_out.dummy_atoms == stats_in.dummy_atoms

    findwaters = FindWaters(refmac.structure, refmac.fphi_best, dummy=True)
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters == stats_in.waters
    assert stats_out.dummy_atoms > stats_in.dummy_atoms

    remove_logs()
