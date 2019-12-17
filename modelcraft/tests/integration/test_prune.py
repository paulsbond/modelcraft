from modelcraft.arguments import parse
from modelcraft.prune import Prune
from modelcraft.gemmineer import model_stats
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
    refmac = Refmac(args, "refmac", data_path("1kv9_model.pdb"))
    prune = Prune("prune", refmac.xyzout, refmac.hklout)
    assert os.path.exists(prune.stdout)
    assert os.path.exists(prune.stderr)
    assert os.path.exists(prune.xyzout)
    statsin = model_stats(refmac.xyzout)
    statsout = model_stats(prune.xyzout)
    assert statsout["residues_built"] > 0
    assert statsout["residues_built"] < statsin["residues_built"]
    os.chdir("..")
    shutil.rmtree(tmp_dir)
