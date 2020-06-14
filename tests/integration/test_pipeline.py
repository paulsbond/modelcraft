import json
import os
import pytest
import shutil
import uuid
from modelcraft.__main__ import main
from modelcraft.reflections import write_mtz
from tests.integration import (
    insulin_contents,
    insulin_fsigf,
    insulin_freer,
    insulin_refmac,
)


def test_insulin():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    fsigf = insulin_fsigf()
    freer = insulin_freer()
    refmac = insulin_refmac()
    write_mtz("data.mtz", [fsigf, freer, refmac.abcd])
    contents = insulin_contents()
    contents.write_sequence_file("sequence.fasta")
    args = []
    args += ["--hklin", "data.mtz"]
    args += ["--seqin", "sequence.fasta"]
    args += ["--cycles", "2"]
    with pytest.raises(SystemExit):
        main(args)
    with open("modelcraft.json") as report_file:
        report = json.load(report_file)
    assert report["real_time"]["total"] > 0
    assert report["termination_reason"] == "Normal"
    os.chdir("..")
    shutil.rmtree(tmp_dir)
