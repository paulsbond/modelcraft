import distutils.spawn
import subprocess
import sys


def run(executable, arguments=[], stdin=[], stdout=None, stderr=None):
    if distutils.spawn.find_executable(executable) is None:
        sys.exit("Executable '%s' not found." % executable)
    p = subprocess.Popen(
        args=[executable] + arguments,
        stdin=subprocess.PIPE if len(stdin) > 0 else None,
        stdout=None if stdout is None else open(stdout, "w"),
        stderr=None if stderr is None else open(stderr, "w"),
        encoding="utf8")
    if len(stdin) > 0:
        for line in stdin:
            p.stdin.write(line + "\n")
        p.stdin.close()
    p.wait()
