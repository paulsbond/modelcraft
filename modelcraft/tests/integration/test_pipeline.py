from modelcraft.tests import data_path
from modelcraft.__main__ import main
import os
import shutil
import uuid


def test_1kv9():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    argument_list = [
        "--hklin", data_path("1kv9_data.mtz"),
        "--seqin", data_path("1kv9_sequence.fasta"),
        "--cycles", "2"
    ]
    main(argument_list)
    os.chdir("..")
    shutil.rmtree(tmp_dir)
