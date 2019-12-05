import json
import time

start_time = time.time()

report = {
    "real_time": {"total": 0},
    "cycles": [],
}


def add_job(job):
    name = job.name.lower().replace(" ", "_")
    if name not in report["real_time"]:
        report["real_time"][name] = 0
    report["real_time"][name] += job.real_time
    write()


def add_cycle(cycle):
    report["cycles"].append(cycle)
    write()


def write():
    report["real_time"]["total"] = round(time.time() - start_time)
    with open("report.json", "w") as report_file:
        json.dump(report, report_file, indent=2)
