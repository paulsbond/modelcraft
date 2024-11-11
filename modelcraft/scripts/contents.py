import argparse
import sys
from ..contents import AsuContents
from ..environ import setup_environ


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
    description = "Create a contents JSON file for a PDB entry"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("entry_id", help="PDB entry ID")
    parser.add_argument("contents", help="Path to write the contents JSON")
    args = parser.parse_args(argument_list)
    contents = AsuContents.from_pdbe(args.entry_id)
    contents.write_json_file(args.contents)


if __name__ == "__main__":
    main()
