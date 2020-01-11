from modelcraft.arguments import parse
from modelcraft.parrot import Parrot
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
    parrot = Parrot(args, "00.00_parrot", args.hklin, args.xyzin)
    assert os.path.exists(parrot.stdout)
    assert os.path.exists(parrot.stderr)
    assert os.path.exists(parrot.hklout.path)
    os.chdir("..")
    shutil.rmtree(tmp_dir)
