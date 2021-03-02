import argparse
import sys
from ..contents import AsuContents


def main(argument_list=None):
    if argument_list is None:
        argument_list = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid")
    parser.add_argument("contents")
    args = parser.parse_args(argument_list)
    contents = AsuContents(args.pdbid)
    contents.write_json_file(args.contents)


if __name__ == "__main__":
    main()
