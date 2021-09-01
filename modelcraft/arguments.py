from typing import List, Optional
import argparse
import os
import gemmi
import numpy
import pandas
from modelcraft.contents import AsuContents
from modelcraft.jobs.freerflag import FreeRFlag
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure


_PARSER = argparse.ArgumentParser(add_help=False)

_req = _PARSER.add_argument_group("Required arguments")
_req.add_argument("--contents", metavar="FILE", required=True)

_opt = _PARSER.add_argument_group("Optional arguments")
_opt.add_argument("--basic", action="store_true")
_opt.add_argument("--convergence-cycles", metavar="N", default=4, type=int)
_opt.add_argument("--convergence-tolerance", metavar="X", default=0.1, type=float)
_opt.add_argument("--cycles", metavar="N", default=25, type=int)
_opt.add_argument("--directory", metavar="PATH", default=".")
_opt.add_argument("--em", action="store_true")
_opt.add_argument("--freerflag", metavar="COL")
_opt.add_argument("--help", action="help")
_opt.add_argument("--keep-jobs", action="store_true")
_opt.add_argument("--keep-logs", action="store_true")
_opt.add_argument("--map", metavar="FILE")
_opt.add_argument("--model", metavar="FILE")
_opt.add_argument("--no-auto-stop", dest="auto_stop", action="store_false")
_opt.add_argument("--observations", metavar="COLS")
_opt.add_argument("--phases", metavar="COLS")
_opt.add_argument("--reflections", metavar="FILE")
_opt.add_argument("--resolution", metavar="X", type=float)
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
    args.contents = AsuContents.from_file(args.contents)
    if args.reflections:
        _parse_data_items(args)
    if args.map:
        _parse_map(args)
        args.basic = True
    if args.model:
        args.model = read_structure(args.model)
    return args


def _basic_check(args: argparse.Namespace):
    if args.em and args.map is None:
        _PARSER.error("A map is required for EM data")
    if args.reflections is None and args.map is None:
        _PARSER.error("Either reflections or a map must be provided")
    if args.reflections is not None and args.map is not None:
        _PARSER.error("Either reflections or a map must be provided (not both)")
    if args.map is None and args.resolution is not None:
        _PARSER.error("The --resolution argument is only used when starting from a map")
    if args.map is not None and args.resolution is None:
        _PARSER.error("The --resolution argument is required when starting from a map")
    if args.resolution is not None and args.resolution <= 0:
        _PARSER.error("Resolution must be greater than 0")
    if args.cycles < 1:
        _PARSER.error("The maximum number of cycles must be greater than 0")
    if args.convergence_cycles < 1:
        _PARSER.error("The number of convergence cycles must be greater than 0")
    if args.convergence_tolerance < 0.1:
        _PARSER.error("The convergence tolerance must be 0.1 or higher")


def _check_paths(args: argparse.Namespace):
    for arg in (
        "buccaneer",
        "contents",
        "map",
        "model",
        "parrot",
        "reflections",
        "sheetbend",
    ):
        path = getattr(args, arg)
        if path is not None:
            path = os.path.abspath(path)
            setattr(args, arg, path)
            if not os.path.exists(path):
                _PARSER.error("File not found: %s" % path)


def _parse_data_items(args: argparse.Namespace):
    mtz = gemmi.read_mtz_file(args.reflections)
    args.observations = _parse_data_item(
        mtz, args.observations, ["FQ", "GLGL", "JQ", "KMKM"], "observations"
    )
    args.freer = _parse_data_item(mtz, args.freerflag, ["I"], "freerflag")
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
        if len(options) == 0:
            _PARSER.error(f"No suitable columns found for the {name}")
        if len(options) > 1:
            message = f"Multiple possible columns found for the {name}:"
            message += "\nPlease add one of the following options:"
            for option in options:
                message += f"\n--{name} {option.label()}"
            _PARSER.error(message)
        print(f"Using {options[0].label()} for the {name}")
        return options[0]
    item = DataItem(mtz, label)
    if item.types not in accepted_types:
        message = f"Column label '{label}' does not match type "
        message += " or ".join(accepted_types)
        _PARSER.error(message)
    return item


def _parse_map(args: argparse.Namespace):
    args.map = gemmi.read_ccp4_map(args.map)
    args.map.setup()
    array = numpy.array(args.map.grid, copy=False)
    if numpy.isnan(array).any():
        _PARSER.error("Map does not cover the full ASU")
    grid = gemmi.transform_map_to_f_phi(args.map.grid, half_l=True)
    data = grid.prepare_asu_data(dmin=args.resolution)
    mtz = gemmi.Mtz()
    mtz.cell = args.map.grid.unit_cell
    mtz.spacegroup = args.map.grid.spacegroup
    mtz.add_dataset("HKL_base")
    mtz.add_column("H", "H")
    mtz.add_column("K", "H")
    mtz.add_column("L", "H")
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    mtz.set_data(data)
    mtz.update_reso()
    array = numpy.array(mtz, copy=True)
    data_frame = pandas.DataFrame(data=array, columns=mtz.column_labels())
    mtz.add_column("SIGF", "Q")
    data_frame["SIGF"] = 1.0
    mtz.add_column("FOM", "W")
    data_frame["FOM"] = 1.0
    mtz.set_data(data_frame.to_numpy())
    args.observations = DataItem(mtz, "F,SIGF")
    args.freer = FreeRFlag(mtz).run().freer
    args.phases = DataItem(mtz, "PHI,FOM")
