from typing import List, Optional
import argparse
import os
import gemmi
from modelcraft.contents import AsuContents
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


_PARSER = argparse.ArgumentParser(add_help=False)

_ = _PARSER.add_argument_group("Required arguments")
_.add_argument("--hklin", metavar="FILE", required=True)
_.add_argument("--seqin", metavar="FILE", required=True)

_ = _PARSER.add_argument_group("Optional arguments")
_.add_argument("--amplitudes", metavar="COLS")
_.add_argument("--convergence-cycles", metavar="N", default=4, type=int)
_.add_argument("--convergence-tolerance", metavar="X", default=0.1, type=float)
_.add_argument("--cycles", metavar="N", default=25, type=int)
_.add_argument("--freerflag", metavar="COL")
_.add_argument("--help", action="help")
_.add_argument("--keep-jobs", action="store_true")
_.add_argument("--model", metavar="FILE")
_.add_argument("--no-auto-stop", dest="auto_stop", action="store_false")
_.add_argument("--remove-non-protein", action="store_true")
_.add_argument("--phases", metavar="COLS")
_.add_argument("--semet", action="store_true")
_.add_argument("--twinned", action="store_true")
_.add_argument("--unbiased", action="store_true")

_ = _PARSER.add_argument_group("Developer arguments")
_.add_argument("--buccaneer", metavar="FILE", default="cbuccaneer")


def parse(arguments: Optional[List[str]] = None) -> argparse.Namespace:
    args = _PARSER.parse_args(arguments)
    _basic_check(args)
    _parse_data_items(args)
    args.contents = AsuContents(args.seqin)
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

    for arg in "hklin", "seqin", "model":
        path = getattr(args, arg)
        if path is not None and not os.path.exists(path):
            _PARSER.error("File not found: %s" % path)


def _parse_data_items(args: argparse.Namespace):
    mtz = gemmi.read_mtz_file(args.hklin)
    args.fsigf = _parse_data_item(mtz, args.amplitudes, ["FQ"], "amplitudes")
    args.freer = _parse_data_item(mtz, args.freerflag, ["I"], "free-R flag")
    if args.phases is not None or args.model is None:
        args.phases = _parse_data_item(mtz, args.phases, ["PW", "AAAA"], "phases")


def _parse_data_item(
    mtz: gemmi.Mtz, label: Optional[str], accepted_types: List[str], name: str
) -> DataItem:
    if label is None:
        options = []
        for types in accepted_types:
            options.extend(DataItem.search(mtz, types))
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
