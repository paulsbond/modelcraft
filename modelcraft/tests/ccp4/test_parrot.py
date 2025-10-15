import gemmi

from modelcraft.contents import AsuContents
from modelcraft.jobs.parrot import Parrot
from modelcraft.monlib import MonLib
from modelcraft.reflections import DataItem

from . import ccp4_path


def test_parrot_after_parrot():
    mtz_path = ccp4_path("examples", "data", "gere.mtz")
    seq_path = ccp4_path("examples", "data", "gere.pir")
    contents = AsuContents.from_file(seq_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    fsigf = DataItem(mtz, "FPHASED,SIGFPHASED")
    freer = DataItem(mtz, "FreeR_flag")
    phases = DataItem(mtz, "HLA,HLB,HLC,HLD")
    monlib = MonLib.STANDARD
    parrot1 = Parrot(contents, fsigf, freer, phases, monlib).run()
    parrot2 = Parrot(contents, fsigf, freer, parrot1.abcd, monlib, parrot1.fphi).run()
    assert parrot2.abcd.nreflections == freer.nreflections
    assert parrot2.fphi.nreflections == freer.nreflections
