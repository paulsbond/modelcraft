import collections
import os
import itertools


class Pipeline:
    def __init__(
        self, directory: str = "", keep_jobs: bool = False, keep_logs: bool = False
    ):
        self._numbers = itertools.count(start=1)
        self.directory = directory
        self.keep_jobs = keep_jobs
        self.keep_logs = keep_logs
        self.seconds = collections.defaultdict(float)

    def path(self, *paths: str) -> str:
        return os.path.join(self.directory, *paths)

    def next_job_directory(self, name: str) -> str:
        return self.path(f"job_{next(self._numbers)}_{name}")
