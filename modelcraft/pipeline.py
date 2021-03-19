import collections
import itertools
from .utils import random_id


class Pipeline:
    def __init__(self, keep_jobs: bool = False, keep_logs: bool = False):
        self._id = random_id(length=10)
        self._numbers = itertools.count(start=1)
        self.keep_jobs = keep_jobs
        self.keep_logs = keep_logs
        self.seconds = collections.defaultdict(float)

    def next_job_directory(self, executable: str) -> str:
        return f"job_{self._id}_{next(self._numbers)}_{executable}"
