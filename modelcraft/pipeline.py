import collections
import itertools


class Pipeline:
    def __init__(self, keep_jobs: bool = False, keep_logs: bool = False):
        self._numbers = itertools.count(start=1)
        self.keep_jobs = keep_jobs
        self.keep_logs = keep_logs
        self.seconds = collections.defaultdict(float)

    def next_job_directory(self, name: str) -> str:
        return f"job_{next(self._numbers)}_{name}"
