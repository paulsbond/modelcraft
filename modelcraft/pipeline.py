import collections
import dataclasses
import itertools
import json
import os
import sys
import time


class Pipeline:
    def __init__(
        self,
        directory: str = "",
        keep_jobs: bool = False,
        keep_logs: bool = False,
        json_name: str = None,
    ):
        self._numbers = itertools.count(start=1)
        self.directory = directory
        self.keep_jobs = keep_jobs
        self.keep_logs = keep_logs
        self.json_name = json_name
        self.seconds = collections.defaultdict(float)
        self.report = {"seconds": self.seconds, "jobs": []}
        self.start_time = None

    def path(self, *paths: str) -> str:
        return os.path.join(self.directory, *paths)

    def next_job_directory(self, name: str) -> str:
        return self.path(f"job_{next(self._numbers)}_{name}")

    def report_job_start(self, name):
        if self.start_time is None:
            self.start_time = time.time()
        print(name, flush=True)
        self.report["running_job"] = name
        self.write_report()

    def report_job_finish(self, result):
        name = self.report.pop("running_job")
        result_dict = {}
        for field in dataclasses.fields(result):
            value = getattr(result, field.name)
            try:
                json.dumps(value)
            except TypeError:
                pass
            else:
                result_dict[field.name] = value
        print(json.dumps(result_dict, indent=4), flush=True)
        self.report["jobs"].append({"name": name, **result_dict})
        self.write_report()

    def write_report(self):
        if self.json_name:
            self.seconds["total"] = time.time() - self.start_time
            with open(self.path(self.json_name), "w") as report_file:
                json.dump(self.report, report_file, indent=4)

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---", flush=True)
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()
