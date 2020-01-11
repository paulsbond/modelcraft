from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
from modelcraft.prune import Prune
from modelcraft.refmac import Refmac
from modelcraft.tests import data_path
import os
import shutil
import uuid


def test_waters():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    argument_list = [
        "--hklin",
        data_path("1kv9_data.mtz"),
        "--seqin",
        data_path("1kv9_sequence.fasta"),
    ]
    args = parse(argument_list)
    xyzin = CoordinateFile(data_path("1kv9_model.pdb"))
    refmac = Refmac(args, "00.01_refmac", xyzin, cycles=1)
    prune = Prune("00.02_prune", refmac.xyzout, refmac.hklout)
    assert os.path.exists(prune.stdout)
    assert os.path.exists(prune.stderr)
    assert os.path.exists(prune.xyzout.path)
    assert prune.xyzout.residues > 0
    assert prune.xyzout.residues < refmac.xyzout.residues
    os.chdir("..")
    shutil.rmtree(tmp_dir)
