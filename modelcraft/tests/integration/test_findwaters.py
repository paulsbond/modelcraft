from modelcraft.arguments import parse
from modelcraft.coordinates import CoordinateFile
from modelcraft.findwaters import FindWaters
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
    findwaters = FindWaters("00.02_findwaters", refmac.xyzout, refmac.hklout)
    assert os.path.exists(findwaters.stdout)
    assert os.path.exists(findwaters.stderr)
    assert os.path.exists(findwaters.xyzout.path)
    assert findwaters.xyzout.waters > 0
    assert findwaters.xyzout.dummys == 0
    finddummys = FindWaters(
        "00.03_finddummys", refmac.xyzout, refmac.hklout, dummy=True
    )
    assert os.path.exists(finddummys.stdout)
    assert os.path.exists(finddummys.stderr)
    assert os.path.exists(finddummys.xyzout.path)
    assert finddummys.xyzout.waters == 0
    assert finddummys.xyzout.dummys > 0
    os.chdir("..")
    shutil.rmtree(tmp_dir)
