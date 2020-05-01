import gemmi
from modelcraft.buccaneer import Buccaneer
from modelcraft.contents import AsuContents
from modelcraft.reflections import FsigF, FreeRFlag, ABCD
from modelcraft.structure import model_stats
from modelcraft.tests import data_path


def test_1kv9():
    contents = AsuContents(data_path("1kv9_sequence.fasta"))
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = FsigF(mtz, "FP,SIGFP")
    freer = FreeRFlag(mtz, "FREE")
    phases = ABCD(mtz, "HL")
    buccaneer = Buccaneer(
        contents=contents, fsigf=fsigf, freer=freer, phases=phases, cycles=1
    )
    stats = model_stats(buccaneer.structure)
    assert stats.residues > 0
