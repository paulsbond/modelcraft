from modelcraft.reflections import DataFile
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
        columns = [
            (args.hklin.path, args.colin_fsigf),
            (args.hklin.path, args.colin_free),
        ]
        if hklin.abcd is not None:
            columns.append((hklin.path, hklin.abcd))
        if hklin.phifom is not None:
            columns.append((hklin.path, hklin.phifom))
        if hklin.fwphiw is not None:
            columns.append((hklin.path, hklin.fwphiw))
        if hklin.fcphic is not None:
            columns.append((hklin.path, hklin.fcphic))
        cmtzjoin_args = ["-mtzout", self.path("hklin.mtz")]
        for path, label in columns:
            cmtzjoin_args += ["-mtzin", path, "-colin", label, "-colout", label]
        self.run("cmtzjoin", cmtzjoin_args)
        new_hklin = DataFile(self.path("hklin.mtz"))
        new_hklin.fsigf = args.colin_fsigf
        new_hklin.free = args.colin_free
        new_hklin.abcd = hklin.abcd
        new_hklin.phifom = hklin.phifom
        new_hklin.fwphiw = hklin.fwphiw
        new_hklin.fcphic = hklin.fcphic
        return new_hklin
