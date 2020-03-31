from modelcraft.__main__ import main
from modelcraft.tests import data_path
import json
import os
import pytest
import shutil
import uuid


def _test_pipeline(argument_list):
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    with pytest.raises(SystemExit):
        main(argument_list)
    with open("modelcraft.json") as f:
        report = json.load(f)
    assert report["real_time"]["total"] > 0
    os.chdir("..")
    shutil.rmtree(tmp_dir)


def test_1kv9_from_phases():
    _test_pipeline(
        ["--hklin", data_path("1kv9_data.mtz"), "--seqin", data_path("1kv9_sequence.fasta"), "--cycles", "2"]
    )


def test_1kv9_from_mr_model():
    _test_pipeline(
        [
            "--hklin",
            data_path("1kv9_data.mtz"),
            "--seqin",
            data_path("1kv9_sequence.fasta"),
            "--mr-model",
            data_path("1kv9_model.pdb"),
            "--cycles",
            "1",
        ]
    )
