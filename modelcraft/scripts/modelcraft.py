import sys
from ..arguments import parse
from ..environ import setup_environ
from ..modelcraftem import ModelCraftEm
from ..modelcraftxray import ModelCraftXray


def main(args=None):
    setup_environ()
    raw_args = args or sys.argv[1:]
    parsed_args = parse(args)
    pipeline = ModelCraftEm if parsed_args.mode == "em" else ModelCraftXray
    pipeline(parsed_args, raw_args).run()


if __name__ == "__main__":
    main()
