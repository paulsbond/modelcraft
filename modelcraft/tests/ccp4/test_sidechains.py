from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi

from modelcraft.reflections import write_mtz
from modelcraft.scripts.sidechains import any_missing_side_chains
from modelcraft.scripts.sidechains import main as fix_side_chains
from modelcraft.structure import write_mmcif

from . import insulin_refmac


def test_insulin():
    refmac = insulin_refmac()
    assert any_missing_side_chains(refmac.structure)
    with TemporaryDirectory() as tempdir:
        xyzin = str(Path(tempdir, "input.cif"))
        hklin = str(Path(tempdir, "input.mtz"))
        xyzout = str(Path(tempdir, "output.cif"))
        write_mmcif(xyzin, refmac.structure)
        write_mtz(hklin, [refmac.fphi_best], ["FWT,PHWT"])
        fix_side_chains([xyzin, hklin, xyzout])
        structure = gemmi.read_structure(xyzout)
        assert not any_missing_side_chains(structure)
