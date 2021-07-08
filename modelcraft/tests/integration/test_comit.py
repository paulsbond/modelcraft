from modelcraft.jobs.comit import Comit
from . import insulin_fsigf, insulin_refmac


def test_insulin():
    fsigf = insulin_fsigf()
    refmac = insulin_refmac()
    comit = Comit(fsigf, refmac.fphi_best).run()
    assert comit.abcd.nreflections == fsigf.nreflections
    assert comit.fphi.nreflections == fsigf.nreflections
