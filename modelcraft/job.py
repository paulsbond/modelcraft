import distutils.spawn
import os
import subprocess
import sys


class Job:
    def __init__(self, directory):
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
            encoding="utf8")
        if len(stdin) > 0:
            for line in stdin:
                p.stdin.write(line + "\n")
            p.stdin.close()
        p.wait()
