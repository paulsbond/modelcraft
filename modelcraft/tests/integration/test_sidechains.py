from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
from modelcraft.refmac import Refmac
from modelcraft.sidechains import Sidechains
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
    refmac1 = Refmac(args, "00.01_refmac", xyzin, cycles=1)
    sidechains = Sidechains("00.02_sidechains", refmac1.xyzout, refmac1.hklout)
    assert os.path.exists(sidechains.stdout)
    assert os.path.exists(sidechains.stderr)
    assert os.path.exists(sidechains.xyzout.path)
    refmac2 = Refmac(args, "00.03_refmac", sidechains.xyzout, cycles=1)
    assert refmac2.final_rfree < refmac1.final_rfree
    os.chdir("..")
    shutil.rmtree(tmp_dir)
