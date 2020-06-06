import os
import sys
from modelcraft.pipeline import Pipeline


def main(argument_list=None):
    if "CCP4" not in os.environ:
        raise EnvironmentError("CCP4 environment not set")
    if argument_list is None:
        argument_list = sys.argv[1:]
    Pipeline(argument_list)


if __name__ == "__main__":
    main()
