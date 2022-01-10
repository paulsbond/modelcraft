from typing import List, Optional
import argparse
import os
import sys
import gemmi
import numpy
import pandas
from . import __version__
from .contents import AsuContents
from .reflections import DataItem
from .structure import read_structure


_PROG = None
if os.path.basename(sys.argv[0]) == "__main__.py":
    _PROG = f"{sys.executable} -m modelcraft"
_PARSER = argparse.ArgumentParser(prog=_PROG)
_PARSER.add_argument("-v", "--version", action="version", version=__version__)

_PARENT = argparse.ArgumentParser(add_help=False)
_GROUP = _PARENT.add_argument_group("required arguments (shared)")
_GROUP.add_argument(
    "--contents",
    required=True,
    metavar="X",
    help=(
        "A file with a description of the assymetric unit contents, "
        "either as a sequence file (in FASTA or PIR format) "
        "with both protein and nucleic acid sequences, "
        "or a more detailed contents file in JSON format. "
        "Example JSON files for existing PDB entries "
        "can be created using the modelcraft-contents script."
    ),
)
_GROUP = _PARENT.add_argument_group("optional arguments (shared)")
_GROUP.add_argument(
    "--model",
    metavar="X",
    help=(
        "A starting model in PDB, mmCIF or mmJSON format. "
        "This could be a placed molecular-replacement model, "
        "a heavy-atom substructure from experimental phasing, "
        "or a partially-built model from another source. "
        "Residues that are not protein, RNA, DNA or water will be kept in place, "
        "so incorporated heavy atoms (such as SE in selenomethionine) "
        "should not be included. "
        "In X-ray mode, the model will first be refined using Sheetbend then Refmac. "
        "If starting phases are not specified using the --phases argument "
        "then phases from the refined model will be used."
    ),
)
_GROUP.add_argument(
    "--cycles",
    default=25,
    type=int,
    metavar="X",
    help="The maximum number of cycles.",
)
_GROUP.add_argument(
    "--auto-stop-cycles",
    default=4,
    type=int,
    metavar="X",
    help=(
        "The number of cycles without improvement "
        "before the program stops automatically. "
        "Improvement is measured by R-free in X-ray mode and FSC in EM mode. "
        "A cycle must improve on the previous best value "
        "to be marked as an improvement. "
        "Setting this value to less than 1 "
        "makes the program run to the maximum number of cycles."
    ),
)
_GROUP.add_argument(
    "--directory",
    default="modelcraft",
    metavar="X",
    help=(
        "The directory where files will be written. "
        "It will be created (along with any intermediate directories) "
        "and can not already exist."
    ),
)
_GROUP.add_argument(
    "--keep-files",
    action="store_true",
    help=(
        "Keep all files from intermediate programs instead of deleting them. "
        "This can lead to large directory sizes."
    ),
)
_GROUP.add_argument(
    "--keep-logs",
    action="store_true",
    help=(
        "Keep log files from intermediate programs "
        "if not using the full --keep-files argument."
    ),
)

_SUB_PARSERS = _PARSER.add_subparsers(dest="mode", required=True)
_FORMATTER = argparse.ArgumentDefaultsHelpFormatter

_XRAY = _SUB_PARSERS.add_parser("xray", parents=[_PARENT], formatter_class=_FORMATTER)
_GROUP = _XRAY.add_argument_group("required arguments (xray)")
_GROUP.add_argument(
    "--data",
    required=True,
    metavar="X",
    help=(
        "Reflection data in MTZ format. "
        "The file must contain merged observations and a free-R flag. "
        "Observations may be amplitudes or intensities, either mean or anomalous. "
        "If a model has not been provided the file must also contain starting phases "
        "as a phase and figure-of-merit or Hendrickson-Lattman coefficients. "
        "An attempt will be made to automatically identify the relevant columns "
        "unless they are specified "
        "using the --observations, --phases and --freerflag arguments."
    ),
)
_GROUP = _XRAY.add_argument_group("optional arguments (xray)")
_GROUP.add_argument(
    "--observations",
    metavar="X",
    dest="observations_label",
    help=(
        "Comma-separated column labels for the observations (e.g. FP,SIGFP). "
        "If anomalous amplitudes or intensities are provided "
        "then CTruncate will be used to convert them to mean amplitudes. "
        "This is not required if the MTZ only contains one set of observations."
    ),
)
_GROUP.add_argument(
    "--phases",
    metavar="X",
    dest="phases_label",
    help=(
        "Comma-separated column labels for the starting phases "
        "as either a phase and figure of merit (e.g. PHIB,FOM) "
        "or Hendrickson-Lattman coefficients (e.g. HLA,HLB,HLC,HLD). "
        "This is not required if the MTZ only contains one set of phases, "
        "unless a model is also provided and the starting phases "
        "should not come from refinement of that model. "
        "Parrot will be used for density modification "
        "before the phases are used for model building."
    ),
)
_GROUP.add_argument(
    "--freerflag",
    metavar="X",
    dest="freerflag_label",
    help=(
        "Column label for the free-R flag (e.g. FreeR_flag). "
        "This is not required if the MTZ only contains one free-R flag."
    ),
)
_GROUP.add_argument(
    "--unbiased",
    action="store_true",
    help=(
        "Pass the starting phases to Refmac for use in refinement "
        "until R-work drops to 35%% or below. "
        "This is usually done after experimental phasing, "
        "but should not be used after molecular-replacement when the phases are biased."
    ),
)
_GROUP.add_argument(
    "--twinned",
    action="store_true",
    help=(
        "Turn on twinned refinement. "
        "Only use this option if you are sure your crystal is twinned."
    ),
)
_GROUP.add_argument(
    "--basic",
    action="store_true",
    help=(
        "Use a more basic pipeline with only Buccaneer, Nautilus and Refmac. "
        "Parrot density modification is still used on the first cycle "
        "and starting models are still refined using Sheetbend and Refmac."
    ),
)

_EM = _SUB_PARSERS.add_parser("em", parents=[_PARENT], formatter_class=_FORMATTER)
_GROUP = _EM.add_argument_group("required arguments (em)")
_GROUP.add_argument(
    "--map",
    required=True,
    metavar="X",
    help="Input map in MRC format",
)
_GROUP.add_argument(
    "--resolution",
    type=float,
    required=True,
    metavar="X",
    help="High resolution limit",
)


def parse(arguments: Optional[List[str]] = None) -> argparse.Namespace:
    args = _PARSER.parse_args(arguments)
    _basic_check(args)
    _check_paths(args)
    args.contents = AsuContents.from_file(args.contents)
    if args.mode == "xray":
        _parse_data_items(args)
    if args.mode == "em":
        _parse_map(args)
    if args.model is not None:
        args.model = read_structure(args.model)
    return args


def _basic_check(args: argparse.Namespace):
    if args.cycles < 1:
        _PARSER.error("--cycles must be greater than 0")
    if args.mode == "em" and args.resolution <= 0:
        _PARSER.error("--resolution must be greater than 0")


def _check_paths(args: argparse.Namespace):
    for arg in ("contents", "data", "map", "model"):
        if hasattr(args, arg):
            path = getattr(args, arg)
            if path is not None:
                if not os.path.exists(path):
                    _PARSER.error("File not found: %s" % path)
                path = os.path.abspath(path)
                setattr(args, arg, path)


def _parse_data_items(args: argparse.Namespace):
    mtz = gemmi.read_mtz_file(args.data)
    _parse_observations(args, mtz)
    _parse_freerflag(args, mtz)
    _parse_phases(args, mtz)


def _parse_observations(args: argparse.Namespace, mtz: gemmi.Mtz):
    label = args.observations_label
    if label is None:
        fmeans = list(DataItem.search(mtz, "FQ"))
        fanoms = list(DataItem.search(mtz, "GLGL"))
        imeans = list(DataItem.search(mtz, "JQ"))
        ianoms = list(DataItem.search(mtz, "KMKM"))
        if any(len(options) > 1 for options in [fmeans, fanoms, imeans, ianoms]):
            _multiple_options_error("observations", fmeans + fanoms + imeans + ianoms)
        if len(fmeans + fanoms + imeans + ianoms) == 0:
            _no_columns_error("observations", ["FQ", "GLGL", "JQ", "KMKM"])
        args.fmean = fmeans[0] if fmeans else None
        args.fanom = fanoms[0] if fanoms else None
        args.imean = imeans[0] if imeans else None
        args.ianom = ianoms[0] if ianoms else None
    else:
        item = _item_from_label(mtz, label, ["FQ", "GLGL", "JQ", "KMKM"])
        args.fmean = item if item.types == "FQ" else None
        args.fanom = item if item.types == "GLGL" else None
        args.imean = item if item.types == "JQ" else None
        args.ianom = item if item.types == "KMKM" else None


def _parse_freerflag(args: argparse.Namespace, mtz: gemmi.Mtz):
    if args.freerflag_label is None:
        freers = list(DataItem.search(mtz, "I"))
        if len(freers) == 0:
            _no_columns_error("free-R flag", ["I"])
        if len(freers) > 1:
            _multiple_options_error("freerflag", freers)
        args.freer = freers[0]
    else:
        args.freer = _item_from_label(mtz, args.freerflag_label, ["I"])
    values = list(args.freer.columns[-1])
    percentage = values.count(0) / len(values) * 100
    if percentage == 0 or percentage > 50:
        _PARSER.error(f"{percentage}% of the reflections are in the free set (flag 0)")


def _parse_phases(args: argparse.Namespace, mtz: gemmi.Mtz):
    if args.phases_label is not None:
        args.phases = _item_from_label(mtz, args.phases_label, ["PW", "AAAA"])
    elif args.model is None:
        phifoms = list(DataItem.search(mtz, "PW", sequential=False))
        abcds = list(DataItem.search(mtz, "AAAA"))
        options = phifoms + abcds
        if len(options) == 0:
            _no_columns_error("phases", ["PW", "AAAA"])
        if len(options) > 1:
            _multiple_options_error("phases", options)
        args.phases = options[0]
    else:
        args.phases = None


def _item_from_label(mtz: gemmi.Mtz, label: str, accepted_types: List[str]) -> DataItem:
    item = DataItem(mtz, label)
    if item.types not in accepted_types:
        message = f"Column label '{item.label()}' does not match type "
        message += " or ".join(accepted_types)
        _PARSER.error(message)
    return item


def _no_columns_error(name: str, types: List[str]):
    message = f"No suitable columns found for the {name}"
    message += f" (with types {' or '.join(types)})"
    _PARSER.error(message)


def _multiple_options_error(argument: str, options: List[DataItem]):
    message = f"Multiple possible columns found for the {argument}:"
    message += "\nPlease add one of the following options:"
    for option in options:
        message += f"\n--{argument} {option.label()}"
    _PARSER.error(message)


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
    args.fmean = DataItem(mtz, "F,SIGF")
    args.phases = DataItem(mtz, "PHI,FOM")
    args.fphi = DataItem(mtz, "F,PHI")
    args.freer = None
