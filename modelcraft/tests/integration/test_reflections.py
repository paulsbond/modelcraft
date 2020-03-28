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
    hklout = refmac.hklout
    assert "FREE,FP,SIGFP,FC,PHIC,FC_ALL,PHIC_ALL,FWT,PHWT,DELFWT,PHDELWT,FOM,FC_ALL_LS,PHIC_ALL_LS" in hklout
    os.chdir("..")
    shutil.rmtree(tmp_dir)
