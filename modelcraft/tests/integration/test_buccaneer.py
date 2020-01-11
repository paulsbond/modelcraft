from modelcraft.arguments import parse
from modelcraft.buccaneer import Buccaneer
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
    ]
    args = parse(argument_list)
    buccaneer = Buccaneer(args, "00.01_buccaneer", args.hklin, cycles=1)
    assert os.path.exists(buccaneer.stdout)
    assert os.path.exists(buccaneer.stderr)
    assert os.path.exists(buccaneer.xyzout.path)
    assert buccaneer.xyzout.residues > 0
    os.chdir("..")
    shutil.rmtree(tmp_dir)
