from modelcraft.reflections import ReflectionFile
import distutils.spawn
import os
import subprocess
import sys
import time


class Job:
    def __init__(self, directory):
        self.start_time = time.time()
        self.cycle = int(directory[:2])
        self.number = int(directory[3:5])
        self.name = directory[6:]
        os.makedirs(directory, exist_ok=True)
        self.directory = os.path.abspath(directory)
        self.stdout = self.path("stdout.txt")
        self.stderr = self.path("stderr.txt")

    def path(self, *paths):
        return os.path.join(self.directory, *paths)

    def run(self, executable, arguments=[], stdin=[]):
        if distutils.spawn.find_executable(executable) is None:
            sys.exit("Executable '%s' not found." % executable)
        p = subprocess.Popen(
            args=[executable] + arguments,
            stdin=subprocess.PIPE if len(stdin) > 0 else None,
            stdout=open(self.stdout, "a"),
            stderr=open(self.stderr, "a"),
            encoding="utf8",
        )
        if len(stdin) > 0:
            for line in stdin:
                p.stdin.write(line + "\n")
            p.stdin.close()
        p.wait()

    def create_hklin(self, args, hklin):
        new_hklin = ReflectionFile(self.path("hklin.mtz"))
        args = [
            "-mtzout",
            new_hklin.path,
            "-mtzin",
            args.hklin.path,
            "-colin",
            args.colin_fsigf,
            "-colout",
            "FP,SIGFP",
            "-mtzin",
            args.hklin.path,
            "-colin",
            args.colin_free,
            "-colout",
            "FREE",
        ]
        new_hklin.fsigf = "FP,SIGFP"
        new_hklin.free = "FREE"
        if hklin.abcd is not None:
            args.extend(
                [
                    "-mtzin",
                    hklin.path,
                    "-colin",
                    hklin.abcd,
                    "-colout",
                    "HLA,HLB,HLC,HLD",
                ]
            )
            new_hklin.abcd = "HLA,HLB,HLC,HLD"
        if hklin.phifom is not None:
            args.extend(
                ["-mtzin", hklin.path, "-colin", hklin.phifom, "-colout", "PHIB,FOM"]
            )
            new_hklin.phifom = "PHIB,FOM"
        if hklin.fphi is not None:
            args.extend(
                ["-mtzin", hklin.path, "-colin", hklin.fphi, "-colout", "FWT,PHWT"]
            )
            new_hklin.fphi = "FWT,PHWT"
        self.run("cmtzjoin", args)
        return new_hklin
