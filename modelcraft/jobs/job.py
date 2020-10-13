from typing import List
import distutils.spawn
import itertools
import os
import shutil
import subprocess
import time


class Job:
    _ids = itertools.count(1)

    def __init__(self, name: str):
        self.id = "%03d_%s" % (next(self._ids), name)
        self._directory = os.path.abspath("job_%s" % self.id)
        self._stdout = self.path("stdout.txt")
        self._stderr = self.path("stderr.txt")
        self._comtxt = self.path("com.txt")

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
            stdout=open(self._stdout, "a"),
            stderr=open(self._stderr, "a"),
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
        with open(self._comtxt, "w") as script_file:
            script_file.write(script)
        os.chmod(self._comtxt, 0o755)

    def finish(self) -> None:
        self.finish_time = time.time()

    def remove(self) -> None:
        logs_dir = "modelcraft-logs"
        os.makedirs(logs_dir, exist_ok=True)
        shutil.copy(self._stdout, os.path.join(logs_dir, "%s.log.txt" % self.id))
        shutil.copy(self._stderr, os.path.join(logs_dir, "%s.err.txt" % self.id))
        shutil.copy(self._comtxt, os.path.join(logs_dir, "%s.com.txt" % self.id))
        shutil.rmtree(self._directory, ignore_errors=True)
