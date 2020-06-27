from modelcraft.jobs import Comit
from tests.integration import remove_logs, insulin_fsigf, insulin_refmac


def test_insulin():
    fsigf = insulin_fsigf()
    refmac = insulin_refmac()
    comit = Comit(fsigf, refmac.fphi_best)
    assert comit.abcd.nreflections == fsigf.nreflections
    assert comit.fphi.nreflections == fsigf.nreflections
    remove_logs()
