from modelcraft.arguments import parse
from modelcraft.comit import Comit
from modelcraft.refmac import Refmac
from modelcraft.tests import data_path
import os
import shutil
import uuid


def test_1kv9():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    argument_list = [
        "--hklin",
        data_path("1kv9_data.mtz"),
        "--seqin",
        data_path("1kv9_sequence.fasta"),
        "--xyzin",
        data_path("1kv9_model.pdb"),
    ]
    args = parse(argument_list)
    refmac = Refmac(args, "00.01_refmac", args.xyzin, cycles=1)
    comit = Comit("00.02_comit", refmac.hklout)
    assert os.path.exists(comit.hklout.path)
    os.chdir("..")
    shutil.rmtree(tmp_dir)
