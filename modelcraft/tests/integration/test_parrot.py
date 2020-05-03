import gemmi
from modelcraft.contents import AsuContents
from modelcraft.jobs import Parrot
from modelcraft.reflections import FsigF, FreeRFlag, ABCD
from modelcraft.tests import data_path


def test_1kv9():
    contents = AsuContents(data_path("1kv9_sequence.fasta"))
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    phases = ABCD(mtz, "HL")
    parrot = Parrot(contents, fsigf, freer, phases)
    assert parrot.abcd.nreflections == mtz.nreflections
    assert parrot.fphi.nreflections == mtz.nreflections
