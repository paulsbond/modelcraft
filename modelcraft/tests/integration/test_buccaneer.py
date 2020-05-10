import gemmi
from modelcraft.jobs import Buccaneer
from modelcraft.contents import AsuContents
from modelcraft.reflections import DataItem
from modelcraft.structure import ModelStats
from modelcraft.tests import data_path


def test_1kv9():
    contents = AsuContents(data_path("1kv9_sequence.fasta"))
    mtz = gemmi.read_mtz_file(data_path("1kv9_data.mtz"))
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    phases = DataItem(mtz, "HL")
    buccaneer = Buccaneer(
        contents=contents, fsigf=fsigf, freer=freer, phases=phases, cycles=1
    )
    stats = ModelStats(buccaneer.structure)
    assert stats.residues > 0
