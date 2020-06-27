from modelcraft.jobs import Buccaneer
from modelcraft.structure import ModelStats
from tests.integration import (
    remove_logs,
    insulin_fsigf,
    insulin_freer,
    insulin_refmac,
    insulin_contents,
)


def test_insulin():
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    refmac = insulin_refmac()
    contents = insulin_contents()
    buccaneer = Buccaneer(contents, fsigf, freer, refmac.abcd, cycles=1)
    stats = ModelStats(buccaneer.structure)
    assert stats.residues > 0
    remove_logs()
