from ..environ import setup_environ
from ..modelcraft import ModelCraft


def main(args=None):
    setup_environ()
    ModelCraft(args).run()


if __name__ == "__main__":
    main()
