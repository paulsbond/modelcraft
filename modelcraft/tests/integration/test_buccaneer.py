from modelcraft.jobs.buccaneer import Buccaneer
from modelcraft.structure import ModelStats
from . import (
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
    buccaneer = Buccaneer(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        freer=freer,
        cycles=1,
    ).run()
    stats = ModelStats(buccaneer.structure)
    assert stats.residues > 0
