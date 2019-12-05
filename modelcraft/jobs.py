import abc
import modelcraft.args as args
import distutils.spawn
import inspect
import itertools
import os
import modelcraft.report as report
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

_initiated_jobs = []


def execute(executable, arguments=[], stdin=[]):
    if distutils.spawn.find_executable(executable) is None:
        sys.exit("Executable '%s' not found." % executable)
    with open("stdout.txt", "a") as stdout, open("stderr.txt", "a") as stderr:
        p = subprocess.Popen(
            args=[executable] + arguments,
            stdin=subprocess.PIPE if len(stdin) > 0 else None,
            stdout=stdout,
            stderr=stderr,
            encoding="utf8")
        if len(stdin) > 0:
            for line in stdin:
                p.stdin.write(line + "\n")
            p.stdin.close()
        p.wait()


def remove_job_directories(cycle):
    if args.keep_intermediate_files:
        return
    for job in _initiated_jobs:
        if job.cycle == cycle:
            shutil.rmtree(job.directory)


class Job(abc.ABC):
    _numbers = itertools.count(1)

    def __init__(self, cycle, name):
        self.cycle = cycle
        self.number = next(self._numbers)
        self.name = name
        self.directory = "job_%d" % (self.number)
        print("%3d %s" % (self.number, name))
        os.mkdir(self.directory)
        self.start_time = time.time()
        _initiated_jobs.append(self)

    def path(self, filename):
        return os.path.join(self.directory, filename)

    def finish(self):
        self.real_time = round(time.time() - self.start_time)
        report.add_job(self)


class Buccaneer(Job):
    def __init__(self, cycle, mtzin, pdbin=None, cycles=2):
        super().__init__(cycle, "Buccaneer")
        self.xmlout = self.path("xmlout.xml")
        self.xyzout = self.path("xyzout.pdb")
        self.run(mtzin, pdbin, cycles)
        self.finish()

    def run(self, mtzin, pdbin, cycles):
        stdin = []
        stdin.append("seqin " + args.seqin)
        stdin.append("colin-fo " + args.colin_fo)
        stdin.append("colin-free " + args.colin_free)
        if mtzin == args.mtzin:
            stdin.append("mtzin " + mtzin)
            if args.colin_hl is not None:
                stdin.append("colin-hl " + args.colin_hl)
            if args.colin_phifom is not None:
                stdin.append("colin-phifom " + args.colin_phifom)
        else:
            execute("cmtzjoin", [
                "-mtzout", self.path("hklin.mtz"),
                "-mtzin", args.mtzin, "-colin", args.colin_fo, "-colout", args.colin_fo,
                "-mtzin", args.mtzin, "-colin", args.colin_free, "-colout", args.colin_free,
                "-mtzin", mtzin, "-colin", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
                                 "-colout", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
                "-mtzin", mtzin, "-colin", "FWT,PHWT", "-colout", "FWT,PHWT",
            ])
            stdin.append("mtzin " + self.path("hklin.mtz"))
            stdin.append("colin-hl HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB")
            stdin.append("colin-fc FWT,PHWT")
        if pdbin is not None:
            stdin.append("pdbin " + pdbin)
            for structure in args.known_structure:
                stdin.append("known-structure " + structure)
        if args.mr_mode > 1:
            stdin.append("pdbin-mr " + args.mr_model)
            if args.mr_mode == 3:
                stdin.append("mr-model")
            if args.mr_mode in (4, 6):
                stdin.append("mr-model-filter")
            if args.mr_mode in (5, 6):
                stdin.append("mr-model-seed")
        stdin.append("pdbout " + self.xyzout)
        stdin.append("xmlout " + self.xmlout)
        stdin.append("anisotropy-correction")
        stdin.append("fast")
        stdin.append("cycles %d" % cycles)
        stdin.append("resolution 2.0")
        stdin.append("model-filter")
        stdin.append("model-filter-sigma 1.0")
        if args.semet:
            stdin.append("build-semet")
        execute(args.buccaneer, ["-stdin"], stdin)
        results = ET.parse(self.path("xmlout.xml")).getroot().find("Final")
        self.completeness_by_residues = float(results.find("CompletenessByResiduesBuilt").text)
        self.completeness_by_chains = float(results.find("CompletenessByChainsBuilt").text)
        self.chains_built = int(results.find("ChainsBuilt").text)
        self.fragments_built = int(results.find("FragmentsBuilt").text)
        self.residues_unique = int(results.find("ResiduesUnique").text)
        self.residues_built = int(results.find("ResiduesBuilt").text)
        self.residues_sequenced = int(results.find("ResiduesSequenced").text)
        self.longest_fragment = int(results.find("ResiduesLongestFragment").text)


class Refmac(Job):
    def __init__(self, cycle, xyzin, use_phases):
        super().__init__(cycle, "REFMAC")
        self.hklout = self.path("hklout.mtz")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, use_phases)
        self.finish()

    def run(self, xyzin, use_phases):
        arguments = [
            "HKLIN", args.mtzin,
            "XYZIN", xyzin,
            "HKLOUT", self.hklout,
            "XYZOUT", self.xyzout,
            "XMLOUT", self.path("xmlout.xml"),
        ]
        stdin = []
        labin = "FP=" + args.colin_fp
        labin += " SIGFP=" + args.colin_sigfp
        labin += " FREE=" + args.colin_free
        if use_phases:
            if args.colin_hl is not None:
                labin += " HLA=" + args.colin_hla
                labin += " HLB=" + args.colin_hlb
                labin += " HLC=" + args.colin_hlc
                labin += " HLD=" + args.colin_hld
            if args.colin_phifom is not None:
                labin += " PHIB=" + args.colin_phi
                labin += " FOM=" + args.colin_fom
        stdin.append("LABIN " + labin)
        stdin.append("NCYCLES 5")
        stdin.append("MAKE HYDR NO")
        if args.twinned:
            stdin.append("TWIN")
        stdin.append("MAKE NEWLIGAND NOEXIT")
        stdin.append("PHOUT")
        stdin.append("PNAME buccaneer")
        stdin.append("DNAME buccaneer")
        stdin.append("END")
        execute("refmac5", arguments, stdin)
        root = ET.parse(self.path("xmlout.xml")).getroot()
        rworks = list(root.iter("r_factor"))
        rfrees = list(root.iter("r_free"))
        self.initial_rwork = float(rworks[0].text)
        self.initial_rfree = float(rfrees[0].text)
        self.final_rwork = float(rworks[-1].text)
        self.final_rfree = float(rfrees[-1].text)


class Prune(Job):
    def __init__(self, cycle, xyzin, hklin, chains_only=False):
        super().__init__(cycle, "Coot pruning")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, hklin, chains_only)
        self.finish()

    def run(self, xyzin, hklin, chains_only):
        coot_prune_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
        coot_prune_path = os.path.join(coot_prune_dir, "coot_prune.py")
        with open(coot_prune_path) as coot_prune:
            script = coot_prune.read()
        script += "deleted = prune(0, 1, 2, %s)\n" % chains_only
        script += "write_pdb_file(0, '%s')\n" % self.xyzout
        script += "import json\n"
        script += "json.dump(deleted, open('%s', 'w'), indent=2)\n" % self.path("deleted.json")
        script += "exit()"
        with open(self.path("script.py"), "w") as script_file:
            script_file.write(script)
        execute("coot", [xyzin, hklin,
                "--no-graphics", "--no-guano", "--script", self.path("script.py")])
        shutil.rmtree("coot-backup", ignore_errors=True)
        shutil.rmtree("coot-download", ignore_errors=True)


class FindWaters(Job):
    def __init__(self, cycle, xyzin, hklin, dummy=False):
        super().__init__(cycle, "Finding waters")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, hklin, dummy)
        self.finish()

    def run(self, xyzin, hklin, dummy):
        arguments = [
            "--pdbin", xyzin,
            "--hklin", hklin,
            "--f", "FWT",
            "--phi", "PHWT",
            "--pdbout", self.path("waters.pdb"),
        ]
        if dummy:
            arguments.append("--flood")
        # --sigma 2.0
        # --min-dist X
        # --max-dist X
        # --flood-atom-radius 1.4 (adjusts contact distance)
        execute("findwaters", arguments)
        execute("pdb_merge", [
            "xyzin1", xyzin,
            "xyzin2", self.path("waters.pdb"),
            "xyzout", self.xyzout,
        ], ["NOMERGE", "END"])
