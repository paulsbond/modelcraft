from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
from modelcraft.sheetbend import Sheetbend
from modelcraft.tests import data_path
import os
import shutil
import uuid


def test_1kv9():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    argument_list = []
    argument_list += ["--hklin", data_path("1kv9_data.mtz")]
    argument_list += ["--seqin", data_path("1kv9_sequence.fasta")]
    args = parse(argument_list)
    xyzin = CoordinateFile(data_path("1kv9_model.pdb"))
    sheetbend = Sheetbend(args, "00.01_sheetbend", xyzin)
    assert os.path.exists(sheetbend.stdout)
    assert os.path.exists(sheetbend.stderr)
    assert os.path.exists(sheetbend.xyzout.path)
    os.chdir("..")
    shutil.rmtree(tmp_dir)
