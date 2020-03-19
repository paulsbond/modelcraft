from modelcraft.coordinates import CoordinateFile
from modelcraft.reflections import (
    ReflectionFile,
    fo_columns,
    free_columns,
    hl_columns,
    phifom_columns,
)
import argparse
import gemmi
import os
import sys


def _argument_parser():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--hklin", metavar="FILE", required=True)
    required.add_argument("--seqin", metavar="FILE", required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--colin-free", metavar="COL")
    optional.add_argument("--colin-fsigf", metavar="COLS")
    optional.add_argument("--colin-hl", metavar="COLS")
    optional.add_argument("--colin-phifom", metavar="COLS")
    optional.add_argument("--cycles", metavar="N", default=25, type=int)
    optional.add_argument("--fix-side-chains", action="store_true")
    optional.add_argument("--free-r-flag", metavar="N", default="0")
    optional.add_argument("--help", action="help")
    optional.add_argument("--keep-intermediate-files", action="store_true")
    optional.add_argument("--known-structure", nargs="+", metavar="SELECTION", default=[])
    optional.add_argument("--mr-mode", metavar="CHOICE", type=int, choices=range(1, 7), default=6)
    optional.add_argument("--mr-model", metavar="FILE")
    optional.add_argument("--no-auto-stop", dest="auto_stop", action="store_false")
    optional.add_argument("--semet", action="store_true")
    optional.add_argument("--twinned", action="store_true")
    optional.add_argument("--unbiased", action="store_true")
    optional.add_argument("--xyzin", metavar="FILE")

    developer = parser.add_argument_group("Developer arguments")
    developer.add_argument("--buccaneer", metavar="FILE", default="cbuccaneer")

    return parser


def _check_paths(args):
    for arg in "hklin", "seqin", "xyzin", "mr_model":
        path = getattr(args, arg)
        if path is None:
            continue
        setattr(args, arg, os.path.abspath(path))
        if not os.path.exists(path):
            print("File not found: %s" % path)
            sys.exit()


def _find_amplitudes(args, mtz):
    options = list(fo_columns(mtz))
    if len(options) == 1:
        print("Using --colin-fsigf %s" % options[0])
        args.colin_fsigf = options[0]
    else:
        if len(options) == 0:
            print("No amplitides found - check input reflection data")
        else:
            print("Multiple amplitudes found - choose one of the following:")
            for option in options:
                print("--colin-fsigf %s" % option)
        sys.exit()


def _find_freer(args, mtz):
    options = list(free_columns(mtz))
    if len(options) == 1:
        print("Using --colin-free %s" % options[0])
        args.colin_free = options[0]
    else:
        if len(options) == 0:
            print("No free-R flag found - check input reflection data.")
        else:
            print("Multiple free-R flags found - choose one of the following:")
            for option in options:
                print("--colin-free %s" % option)
        sys.exit()


def _find_phases(args, mtz):
    hl_options = list(hl_columns(mtz))
    phifom_options = list(phifom_columns(mtz))
    if len(hl_options) + len(phifom_options) == 1:
        if len(hl_options) == 1:
            print("Using --colin-hl %s" % hl_options[0])
            args.colin_hl = hl_options[0]
        else:
            print("Using --colin-phifom %s" % phifom_options[0])
            args.colin_phifom = phifom_options[0]
    else:
        if len(hl_options) + len(phifom_options) == 0:
            print("No phases found - check input reflection data.")
        else:
            print("Multiple phases found - choose one of the following:")
            for hl_option in hl_options:
                print("--colin-hl %s" % hl_option)
            for phifom_option in phifom_options:
                print("--colin-phifom %s" % phifom_option)
        sys.exit()


def _find_mtz_columns(args):
    mtz = gemmi.read_mtz_file(args.hklin)

    if args.colin_fsigf is None:
        print("\nInput amplitudes not provided")
        _find_amplitudes(args, mtz)

    if args.colin_free is None:
        print("\nInput free-R flag not provided")
        _find_freer(args, mtz)

    if args.colin_hl is None and args.colin_phifom is None and args.mr_model is None:
        print("\nInput phases not provided")
        _find_phases(args, mtz)


def _derive_other_args(args):
    args.hklin = ReflectionFile(args.hklin, args.colin_fsigf, args.colin_free, args.colin_hl, args.colin_phifom)
    if args.xyzin is not None:
        args.xyzin = CoordinateFile(args.xyzin)
    if args.mr_model is not None:
        args.mr_model = CoordinateFile(args.mr_model)

    fo = args.colin_fsigf.split(",")
    args.colin_fp = fo[0]
    args.colin_sigfp = fo[1]

    if args.colin_hl is not None:
        hl = args.colin_hl.split(",")
        args.colin_hla = hl[0]
        args.colin_hlb = hl[1]
        args.colin_hlc = hl[2]
        args.colin_hld = hl[3]

    if args.colin_phifom is not None:
        phifom = args.colin_phifom.split(",")
        args.colin_phi = phifom[0]
        args.colin_fom = phifom[1]


def _check(args):
    if args.cycles < 1:
        print("The maximum number of cycles must be greater than 0")
        sys.exit()

    if args.colin_hl is not None and args.colin_phifom is not None:
        print("Two different input phases are provided. Please choose one:")
        print("--colin-hl %s" % args.colin_hl)
        print("--colin-phifom %s" % args.colin_phifom)
        sys.exit()


def parse(argument_list):
    print("\n %s" % " ".join(argument_list).replace(" --", "\n --"))
    parser = _argument_parser()
    args = parser.parse_args(argument_list)
    _check_paths(args)
    _find_mtz_columns(args)
    _derive_other_args(args)
    _check(args)
    return args
