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
        self.id = _generate_id()
        self._directory = os.path.abspath("job_%s" % self.id)
        self._stdout_path = self.path("stdout.txt")
        self._stderr_path = self.path("stderr.txt")
        self._comtxt_path = self.path("com.txt")
        self.stdout: str = ""
        self.stderr: str = ""
        self.comtxt: str = ""

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
        self._write_cmd_script(executable, arguments, stdin)
        process = subprocess.Popen(
            args=[executable] if arguments is None else ([executable] + arguments),
            stdin=None if stdin is None else subprocess.PIPE,
            stdout=open(self._stdout_path, "a"),
            stderr=open(self._stderr_path, "a"),
            encoding="utf8",
        )
        if stdin is not None:
            for line in stdin:
                process.stdin.write(line + "\n")
            process.stdin.close()

        process.wait()

    def _write_cmd_script(
        self, executable: str, arguments: List[str] = None, stdin: List[str] = None
    ) -> None:
        script = "#!/usr/bin/env bash\n\n"
        script += executable
        if arguments is not None:
            for argument in arguments:
                if self._directory + "/" in argument:
                    argument = argument.split(self._directory + "/")[-1]
                script += " " + argument
        script += " \\\n> stdout.txt 2> stderr.txt"
        if stdin is not None:
            script += " << EOF\n"
            for line in stdin:
                script += f"{line}\n"
            script += "EOF\n"
        else:
            script += "\n"
        with open(self._comtxt_path, "w") as script_file:
            script_file.write(script)
        os.chmod(self._comtxt_path, 0o755)

    def finish(self) -> None:
        self.finish_time = time.time()
        with open(self._stdout_path) as f:
            self.stdout = f.read()
        with open(self._stderr_path) as f:
            self.stderr = f.read()
        with open(self._comtxt_path) as f:
            self.comtxt = f.read()
        shutil.rmtree(self._directory, ignore_errors=True)
