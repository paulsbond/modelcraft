from modelcraft.jobs import Parrot
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
    parrot = Parrot(contents, fsigf, freer, refmac.abcd)
    assert parrot.abcd.nreflections == fsigf.nreflections
    assert parrot.fphi.nreflections == fsigf.nreflections
    remove_logs()
