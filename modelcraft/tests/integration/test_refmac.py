from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
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
    ]
    args = parse(argument_list)
    xyzin = CoordinateFile(data_path("1kv9_model.pdb"))
    refmac = Refmac(args, "00.01_refmac", xyzin, cycles=1)
    assert os.path.exists(refmac.stdout)
    assert os.path.exists(refmac.stderr)
    assert os.path.exists(refmac.hklout.path)
    assert os.path.exists(refmac.xmlout)
    assert os.path.exists(refmac.xyzout.path)
    assert refmac.final_rfree < 100
    assert refmac.final_rwork < refmac.initial_rwork
    os.chdir("..")
    shutil.rmtree(tmp_dir)
