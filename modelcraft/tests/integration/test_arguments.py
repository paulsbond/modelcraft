from modelcraft.arguments import parse
from . import ccp4_path


def test1():
    seqin = ccp4_path("examples", "data", "1mzr_1a80.pir")
    hklin = ccp4_path("examples", "data", "1rxf.mtz")
    xyzin = ccp4_path("examples", "data", "1rxf_randomise.pdb")
    args = ["xray"]
    args += ["--contents", seqin]
    args += ["--data", hklin]
    args += ["--model", xyzin]
    parse(args)


def test2():
    seqin = ccp4_path("examples", "data", "gere.seq")
    hklin = ccp4_path("examples", "data", "gere.mtz")
    args = ["xray"]
    args += ["--contents", seqin]
    args += ["--data", hklin]
    args += ["--observations", "FPHASED,SIGFPHASED"]
    args += ["--phases", "HLA,HLB,HLC,HLD"]
    args += ["--freerflag", "FreeR_flag"]
    args += ["--cycles", "100"]
    args += ["--auto-stop-cycles", "2"]
    args += ["--directory", "tmp"]
    args += ["--keep-files"]
    args += ["--unbiased"]
    args += ["--twinned"]
    args += ["--basic"]
    parse(args)
