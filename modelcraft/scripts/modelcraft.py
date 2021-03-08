import sys
from ..environ import setup_environ
from ..pipeline import Pipeline


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
    Pipeline(argument_list)


if __name__ == "__main__":
    main()
