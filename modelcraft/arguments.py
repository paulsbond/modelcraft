from typing import List, Optional
import argparse
import os
import gemmi
from modelcraft.contents import AsuContents
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


_PARSER = argparse.ArgumentParser(add_help=False)

_req = _PARSER.add_argument_group("Required arguments")
_req.add_argument("--data", metavar="FILE", required=True)
_req.add_argument("--contents", metavar="FILE", required=True)

_opt = _PARSER.add_argument_group("Optional arguments")
_opt.add_argument("--basic", action="store_true")
_opt.add_argument("--convergence-cycles", metavar="N", default=4, type=int)
_opt.add_argument("--convergence-tolerance", metavar="X", default=0.1, type=float)
_opt.add_argument("--cycles", metavar="N", default=25, type=int)
_opt.add_argument("--directory", metavar="PATH", default=".")
_opt.add_argument("--freerflag", metavar="COL")
_opt.add_argument("--help", action="help")
_opt.add_argument("--keep-jobs", action="store_true")
_opt.add_argument("--keep-logs", action="store_true")
_opt.add_argument("--model", metavar="FILE")
_opt.add_argument("--no-auto-stop", dest="auto_stop", action="store_false")
_opt.add_argument("--observations", metavar="COLS")
_opt.add_argument("--phases", metavar="COLS")
_opt.add_argument("--twinned", action="store_true")
_opt.add_argument("--unbiased", action="store_true")

_dev = _PARSER.add_argument_group("Developer arguments")
_dev.add_argument("--buccaneer", metavar="FILE")
_dev.add_argument("--parrot", metavar="FILE")
_dev.add_argument("--sheetbend", metavar="FILE")


def parse(arguments: Optional[List[str]] = None) -> argparse.Namespace:
    args = _PARSER.parse_args(arguments)
    _basic_check(args)
    _check_paths(args)
    _parse_data_items(args)
    args.contents = AsuContents.from_file(args.contents)
    if args.model is not None:
        args.model = read_structure(args.model)
    return args


def _basic_check(args: argparse.Namespace):
    if args.cycles < 1:
        _PARSER.error("The maximum number of cycles must be greater than 0")

    if args.convergence_cycles < 1:
        _PARSER.error("The number of convergence cycles must be greater than 0")

    if args.convergence_tolerance < 0.1:
        _PARSER.error("The convergence tolerance must be 0.1 or higher")


def _check_paths(args: argparse.Namespace):
    for arg in "data", "contents", "model", "buccaneer", "parrot", "sheetbend":
        path = getattr(args, arg)
        if path is not None:
            path = os.path.abspath(path)
            setattr(args, arg, path)
            if not os.path.exists(path):
                _PARSER.error("File not found: %s" % path)


def _parse_data_items(args: argparse.Namespace):
    mtz = gemmi.read_mtz_file(args.data)
    args.observations = _parse_data_item(
        mtz, args.observations, ["FQ", "GLGL", "JQ", "KMKM"], "observations"
    )
    args.freer = _parse_data_item(mtz, args.freerflag, ["I"], "free-R flag")
    if args.phases is not None or args.model is None:
        args.phases = _parse_data_item(mtz, args.phases, ["PW", "AAAA"], "phases")


def _parse_data_item(
    mtz: gemmi.Mtz, label: Optional[str], accepted_types: List[str], name: str
) -> DataItem:
    if label is None:
        options = []
        for types in accepted_types:
            items = DataItem.search(mtz, types, sequential=(types != "PW"))
            options.extend(items)
        if len(options) == 1:
            print(f"Using {options[0].label()} for the {name}")
            return options[0]
        if len(options) > 1:
            message = f"Multiple possible columns found for the {name}:"
            for option in options:
                message += f"\n{option.label()}"
            _PARSER.error(message)
        _PARSER.error(f"No suitable columns found for the {name}")
    else:
        item = DataItem(mtz, label)
        if item.types in accepted_types:
            return item
        message = f"Column label '{label}' does not match type "
        message += " or ".join(accepted_types)
        _PARSER.error(message)
