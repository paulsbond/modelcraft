import gemmi
from modelcraft.reflections import DataItem
from modelcraft.contents import AsuContents
from modelcraft.jobs.parrot import Parrot
from . import ccp4_path


def test_parrot_after_parrot():
    mtz_path = ccp4_path("examples", "data", "gere.mtz")
    seq_path = ccp4_path("examples", "data", "gere.pir")
    contents = AsuContents.from_file(seq_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "FPHASED,SIGFPHASED")
    freer = DataItem(mtz, "FreeR_flag")
    phases = DataItem(mtz, "HLA,HLB,HLC,HLD")
    parrot1 = Parrot(contents, fsigf, freer, phases).run()
    parrot2 = Parrot(contents, fsigf, freer, parrot1.abcd, parrot1.fphi).run()
    assert parrot2.abcd.nreflections == freer.nreflections
    assert parrot2.fphi.nreflections == freer.nreflections
