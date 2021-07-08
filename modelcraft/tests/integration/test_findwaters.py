from modelcraft.jobs.findwaters import FindWaters
from modelcraft.structure import ModelStats
from . import insulin_refmac


def test_insulin_water():
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


def test_insulin_dummy():
    refmac = insulin_refmac()
    stats_in = ModelStats(refmac.structure)
    findwaters = FindWaters(
        structure=refmac.structure,
        fphi=refmac.fphi_best,
        dummy=True,
    ).run()
    stats_out = ModelStats(findwaters.structure)
    assert stats_out.residues == stats_in.residues
    assert stats_out.waters == stats_in.waters
    assert stats_out.dummy_atoms > stats_in.dummy_atoms


def test_existing_water_chain():
    refmac = insulin_refmac()
    findwaters1 = FindWaters(
        structure=refmac.structure,
        fphi=refmac.fphi_best,
    ).run()
    findwaters2 = FindWaters(
        structure=findwaters1.structure,
        fphi=refmac.fphi_best,
    ).run()
    chains1 = len(findwaters1.structure[0])
    chains2 = len(findwaters2.structure[0])
    assert chains1 == chains2
    waters1 = ModelStats(findwaters1.structure).waters
    waters2 = ModelStats(findwaters2.structure).waters
    assert waters2 > waters1
