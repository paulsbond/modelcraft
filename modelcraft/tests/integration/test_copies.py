import os
import shutil
import uuid
from modelcraft.scripts.copies import main
from . import ccp4_path, insulin_contents


def test_insulin():
    tmp_dir = "tmp%s" % uuid.uuid4()
    os.mkdir(tmp_dir)
    os.chdir(tmp_dir)
    contents = insulin_contents()
    contents.write_json_file("contents.json")
    args = ["contents.json", ccp4_path("examples", "data", "insulin.mtz")]
    main(args)
    os.chdir("..")
    shutil.rmtree(tmp_dir)
