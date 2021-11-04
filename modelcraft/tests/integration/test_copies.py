from modelcraft.scripts.copies import main
from . import ccp4_path, in_temp_directory, insulin_contents


@in_temp_directory
def test_insulin():
    contents = insulin_contents()
    contents.write_json_file("contents.json")
    args = ["contents.json", ccp4_path("examples", "data", "insulin.mtz")]
    main(args)
