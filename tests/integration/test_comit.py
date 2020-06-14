from modelcraft.jobs import Comit
from tests.integration import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    Comit(refmac.fsigf, refmac.fphi_best)
