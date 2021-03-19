import sys
from ..environ import setup_environ
from ..modelcraft import ModelCraft


def main(argument_list=None):
    setup_environ()
    if argument_list is None:
        argument_list = sys.argv[1:]
    ModelCraft(argument_list).run()


if __name__ == "__main__":
    main()
