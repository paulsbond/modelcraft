import abc
import os
import distutils.spawn
import shutil
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
import gemmi
from .pipeline import Pipeline
from .reflections import write_mtz
from .structure import read_structure, write_mmcif


class Job(abc.ABC):
    def __init__(self, executable: str):
        self._executable = executable
        self._args = []
        self._stdin = []
        self._environ = {}
        self._cifins = {}
        self._hklins = {}
        self._xyzins = {}
        self._cifouts = {}
        self._hklouts = {}
        self._xmlouts = {}
        self._xyzouts = {}
        self._directory = None
        self._seconds = None

    @abc.abstractmethod
    def run(self, pipeline: Pipeline = None) -> None:
        if distutils.spawn.find_executable(self._executable) is None:
            raise ValueError("Executable '%s' not found" % self._executable)
        self._write_files(pipeline)
        self._run_subprocess()
        self._read_files()
        if pipeline is None:
            self._remove_files()
        else:
            pipeline.times.setdefault(self._executable, 0)
            pipeline.times[self._executable] += self._seconds
            if not pipeline.keep_jobs:
                self._remove_files(keep_logs=pipeline.keep_logs)

    def _path(self, *paths: str) -> str:
        return os.path.join(self._directory, *paths)

    def _write_files(self, pipeline: Pipeline = None) -> str:
        if pipeline is None:
            self._directory = str(uuid.uuid4())
        else:
            self._directory = pipeline.directory(self._executable)
        os.mkdir(self._directory)
        with open(self._path("script.sh"), "w") as stream:
            stream.write(self._script())
        os.chmod(self._path("script.sh"), 0o755)
        for filename, items in self._hklins.items():
            write_mtz(self._path(filename), items)
        for filename, structure in self._xyzins.items():
            write_mmcif(self._path(filename), structure)
        for filename, document in self._cifins.items():
            document.write_file(self._path(filename))

    def _run_subprocess(self):
        start_time = time.time()
        with open(self._path("stdout.txt"), "w") as out_stream:
            with open(self._path("stderr.txt"), "w") as err_stream:
                process = subprocess.Popen(
                    args=[self._executable] + self._args,
                    stdin=subprocess.PIPE if self._stdin else None,
                    stdout=out_stream,
                    stderr=err_stream,
                    encoding="utf8",
                    env={**os.environ, **self._environ},
                    cwd=self._directory,
                )
        if self._stdin:
            for line in self._stdin:
                process.stdin.write(line + "\n")
            process.stdin.close()
        process.wait()
        self._seconds = time.time() - start_time

    def _read_files(self):
        for filename in self._cifouts:
            self._cifouts[filename] = gemmi.cif.read(self._path(filename))
        for filename in self._hklouts:
            self._hklouts[filename] = gemmi.read_mtz_file(self._path(filename))
        for filename in self._xmlouts:
            self._xmlouts[filename] = ET.parse(self._path(filename)).getroot()
        for filename in self._xyzouts:
            self._xyzouts[filename] = read_structure(self._path(filename))

    def _script(self) -> str:
        script = "#!/usr/bin/env bash\n\n"
        if self._environ:
            for variable, value in self._environ.items():
                script += f"export {variable}={value}\n"
            script += "\n"
        script += self._executable
        script += f" {' '.join(self._args)} \\\n> stdout.txt 2> stderr.txt"
        if self._stdin:
            script += " << EOF\n"
            for line in self._stdin:
                script += f"{line}\n"
            script += "EOF\n"
        else:
            script += "\n"
        return script

    def _remove_files(self, keep_logs: bool = False) -> None:
        if keep_logs:
            os.makedirs("modelcraft-logs", exist_ok=True)
            for filename in ("stdout.txt", "stderr.txt", "script.sh"):
                src = self._path(filename)
                dst = os.path.join("modelcraft-logs", f"{self._directory}_{filename}")
                os.rename(src, dst)
        shutil.rmtree(self._directory, ignore_errors=True)
