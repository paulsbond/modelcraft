from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
from modelcraft.findwaters import FindWaters
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
        "--hklin", data_path("1kv9_data.mtz"),
        "--seqin", data_path("1kv9_sequence.fasta"),
    ]
    args = parse(argument_list)
    xyzin = CoordinateFile(data_path("1kv9_model.pdb"))
    refmac = Refmac(args, "refmac", xyzin)
    findwaters = FindWaters("findwaters", refmac.xyzout, refmac.hklout)
    assert os.path.exists(findwaters.stdout)
    assert os.path.exists(findwaters.stderr)
    assert os.path.exists(findwaters.xyzout.path)
    assert findwaters.xyzout.waters > 1
    os.chdir("..")
    shutil.rmtree(tmp_dir)
