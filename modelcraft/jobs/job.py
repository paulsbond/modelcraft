from typing import List
import distutils.spawn
import os
import random
import shutil
import string
import subprocess
import time


def _generate_id():
    chars = string.ascii_lowercase + string.digits
    return "".join(random.choice(chars) for _ in range(8))


class Job:
    def __init__(self):
        self._directory = os.path.abspath("job_%s" % _generate_id())
        self._stdout = self.path("stdout.txt")
        self._stderr = self.path("stderr.txt")
        os.mkdir(self._directory)
        self.start_time = time.time()
        self.finish_time = None

    @property
    def name(self) -> str:
        return type(self).__name__

    def path(self, *paths: str) -> str:
        return os.path.join(self._directory, *paths)

    def run(
        self, executable: str, arguments: List[str] = None, stdin: List[str] = None
    ) -> None:
        if distutils.spawn.find_executable(executable) is None:
            raise ValueError("Executable '%s' not found" % executable)
        process = subprocess.Popen(
            args=[executable] if arguments is None else ([executable] + arguments),
            stdin=None if stdin is None else subprocess.PIPE,
            stdout=open(self._stdout, "a"),
            stderr=open(self._stderr, "a"),
            encoding="utf8",
        )
        if stdin is not None:
            for line in stdin:
                process.stdin.write(line + "\n")
            process.stdin.close()
        process.wait()

    def finish(self) -> None:
        self.finish_time = time.time()
        shutil.rmtree(self._directory, ignore_errors=True)
