from modelcraft.jobs import FindWaters
from modelcraft.structure import ModelStats
from tests.integration import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    stats_in = ModelStats(refmac.structure)

    findwaters = FindWaters(
        structure=refmac.structure,
        fphi=refmac.fphi_best,
    ).run()
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters > stats_in.waters
    assert stats_out.dummy_atoms == stats_in.dummy_atoms

    findwaters = FindWaters(
        structure=refmac.structure,
        fphi=refmac.fphi_best,
        dummy=True,
    ).run()
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters == stats_in.waters
    assert stats_out.dummy_atoms > stats_in.dummy_atoms
