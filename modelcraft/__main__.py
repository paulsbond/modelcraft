from modelcraft.pipeline import Pipeline
import sys


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    Pipeline(args)


if __name__ == "__main__":
    main()
