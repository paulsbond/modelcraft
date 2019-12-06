import argparse
import gemmi
import modelcraft.gemmineer as gemmineer
import sys

# Placeholders for linting and to avoid duplicate declarations
add_waters = None
auto_stop = None
buccaneer = None
colin_fo = None
colin_fom = None
colin_fp = None
colin_free = None
colin_hl = None
colin_hla = None
colin_hlb = None
colin_hlc = None
colin_hld = None
colin_phi = None
colin_phifom = None
colin_sigfp = None
cycles = None
free_r_flag = None
keep_intermediate_files = None
known_structure = []
mr_mode = None
mr_model = None
mtzin = None
semet = None
seqin = None
twinned = None
unbiased = None
xyzin = None


def parse(arguments):
    print("\n %s" % " ".join(arguments).replace(" --", "\n --"))
    args = argument_parser().parse_args(arguments)
    find_mtz_columns(args)
    derive_other_args(args)
    check(args)
    globals().update(vars(args))


def argument_parser():
    parser = argparse.ArgumentParser(
        description="ModelCraft: An automated model building pipeline for X-ray crystallography and cryo-EM",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    required = parser.add_argument_group("Required arguments")
    required.add_argument("--mtzin", metavar="FILE", required=True, help="Input reflection data in MTZ format")
    required.add_argument("--seqin", metavar="FILE", required=True, help="Input protein sequence")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--add-waters", action="store_true",
                          help="Add water molecules once R-work drops below 0.4")
    optional.add_argument("--buccaneer", metavar="FILE", default="cbuccaneer",
                          help="Path to an alternative buccaneer binary")
    optional.add_argument("--colin-fo", metavar="COLS",
                          help=("Column labels for the observed amplitudes\n"
                                "(e.g. FP,SIGFP)"))
    optional.add_argument("--colin-free", metavar="COL",
                          help=("Column label for the free-R flag\n"
                                "(e.g. FreeR_flag)"))
    optional.add_argument("--colin-hl", metavar="COLS",
                          help=("Column labels for input phases as Hendrickson-Lattman coefficients\n"
                                "(e.g. HLA,HLB,HLC,HLD)"))
    optional.add_argument("--colin-phifom", metavar="COLS",
                          help=("Column labels for input phases as a phase and figure of merit\n"
                                "(e.g. PHIB,FOM)"))
    optional.add_argument("--cycles", metavar="N", default=25, type=int,
                          help="Maximum number of pipeline cycles")
    optional.add_argument("--free-r-flag", metavar="N", default="0",
                          help="Flag that identifies the free reflections (default: 0)")
    optional.add_argument("--help", action="help",
                          help="Show this help message and exit")
    optional.add_argument("--keep-intermediate-files", action="store_true",
                          help="Don't delete intermediate files")
    optional.add_argument("--known-structure", nargs="+", metavar="SELECTION", default=[],
                          help=("Known structure selections from the input coordinates\n"
                                "Buccaneer will avoid building into these selections\n"
                                "and they will be copied into the output coordinates\n"
                                "Multiple selections can be specified\n"
                                "The format is /[chain]/[residue]/[atom]:[radius]\n"
                                "e.g. /A/*/*/:2.0         avoid building within 2A of the A chain\n"
                                "e.g. /*/*/ZN    /:3.0    avoid building within 3A of ZN atoms\n"))
    optional.add_argument("--mr-mode", metavar="CHOICE", type=int, choices=range(1, 7), default=1,
                          help=("Determine how the molecular replacement model is used:\n"
                                "1 - Phasing\n"
                                "2 - Phasing, placing/naming chains\n"
                                "3 - Phasing, placing/naming chains, copy every residue\n"
                                "4 - Phasing, placing/naming chains, copy every filtered residue\n"
                                "5 - Phasing, placing/naming chains, copy every 3rd residue\n"
                                "6 - Phasing, placing/naming chains, copy every 3rd filtered residue\n"
                                "(default: 1)"))
    optional.add_argument("--mr-model", metavar="FILE",
                          help=("Input placed molecular replacement model\n"
                                "If input phases are not specified this will be refined"))
    optional.add_argument("--no-auto-stop", dest="auto_stop", action="store_false",
                          help="Run the maximum number of cycles even if the model is not improving")
    optional.add_argument("--semet", action="store_true",
                          help="Build selenomethionine instead of methionine")
    optional.add_argument("--twinned", action="store_true",
                          help=("Turn on twinned refinement\n"
                                "Only do this if you are sure your data are twinned"))
    optional.add_argument("--unbiased", action="store_true",
                          help="Pass input phases to REFMAC for MLHL refinement")
    optional.add_argument("--xyzin", metavar="FILE",
                          help="Input starting coordinates")

    return parser


def find_mtz_columns(args):
    mtz = gemmi.read_mtz_file(args.mtzin)

    if args.colin_fo is None:
        print("\nInput amplitudes not provided")
        options = list(gemmineer.fo_columns(mtz))
        if len(options) == 1:
            print("Using --colin-fo %s" % options[0])
            args.colin_fo = options[0]
        else:
            if len(options) == 0:
                print("No amplitides found - check input reflection data")
            else:
                print("Multiple amplitudes found - choose one of the following:")
                for option in options:
                    print("--colin-fo %s" % option)
            sys.exit()

    if args.colin_free is None:
        print("\nInput free-R flag not provided")
        options = list(gemmineer.free_columns(mtz))
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

    if args.colin_hl is None and args.colin_phifom is None and args.mr_model is None:
        print("\nInput phases not provided")
        hl_options = list(gemmineer.hl_columns(mtz))
        phifom_options = list(gemmineer.phifom_columns(mtz))
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


def derive_other_args(args):
    fo = args.colin_fo.split(",")
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


def check(args):
    if args.cycles < 1:
        print("The maximum number of cycles must be greater than 0")
        sys.exit()

    if args.colin_hl is not None and args.colin_phifom is not None:
        print("Two different input phases are provided. Please choose one:")
        print("--colin-hl %s" % args.colin_hl)
        print("--colin-phifom %s" % args.colin_phifom)
        sys.exit()
