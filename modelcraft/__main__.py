from modelcraft.pipeline import Pipeline
import os
import sys


def main(args=None):
    if "CCP4" not in os.environ:
        sys.exit("CCP4 environment not set")
    if args is None:
        args = sys.argv[1:]
    Pipeline(args)


if __name__ == "__main__":
    main()
