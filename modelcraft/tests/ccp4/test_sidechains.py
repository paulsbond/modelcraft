from modelcraft.jobs.sidechains import Sidechains
from modelcraft.scripts.sidechains import any_missing_side_chains

from . import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    assert any_missing_side_chains(refmac.structure)
    sidechains = Sidechains(refmac.structure, refmac.fphi_best).run()
    assert not any_missing_side_chains(sidechains.structure)
