import gemmi
from modelcraft.jobs.phasematch import PhaseMatch
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import ccp4_path


def test_gere_phases():
    mtz = gemmi.read_mtz_file(ccp4_path("examples", "data", "gere.mtz"))
    fsigf = DataItem(mtz, "FPHASED,SIGFPHASED")
    phases1 = DataItem(mtz, "HLA,HLB,HLC,HLD")
    phases2 = DataItem(mtz, "PHIDM,FOMDM")
    phasematch = PhaseMatch(fsigf, phases1, phases2).run()
    assert 0 < phasematch.f_map_correlation < 1


def test_gere_refmac():
    mtz = gemmi.read_mtz_file(ccp4_path("examples", "data", "gere.mtz"))
    fsigf = DataItem(mtz, "FPHASED,SIGFPHASED")
    freer = DataItem(mtz, "FreeR_flag")
    model1 = read_structure(ccp4_path("examples", "data", "gere_heavy.pdb"))
    model2 = read_structure(ccp4_path("examples", "data", "gere_incompl.pdb"))
    refmac1 = RefmacXray(model1, fsigf, freer, cycles=0).run()
    refmac2 = RefmacXray(model2, fsigf, freer, cycles=0).run()
    phasematch = PhaseMatch(fsigf, refmac1.abcd, refmac2.abcd).run()
    assert 0 < phasematch.f_map_correlation < 1
