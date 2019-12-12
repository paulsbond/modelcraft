import abc
import distutils.spawn
import itertools
import json
import modelcraft.gemmineer as gemmineer
import os
import shutil
import subprocess
import sys
import time


def add_cycle(cycle):
    report["cycles"].append(cycle)
    write()


def write():
    report["real_time"]["total"] = round(time.time() - start_time)
    with open("report.json", "w") as report_file:
        json.dump(report, report_file, indent=2)


class Pipeline(abc.ABC):
    def __init__(self, directory=None):
        self.directory = os.path.abspath(os.curdir if directory is None else directory)
        os.makedirs(self.directory, exist_ok=True)
        os.chdir(self.directory)
        self.start_time = time.time()
        self.cycles = []
        self.jobs = []
        self.real_time_taken = {"total": 0}

    class Job(abc.ABC):
        def __init__(self, cycle, name):
            self.dir
            self.cycle = cycle
            self.number = next(self._numbers)
            self.name = name
            self.directory = "job_%d" % (self.number)
            print("%3d %s" % (self.number, name))
            os.mkdir(self.directory)
            self.start_time = time.time()
            _initiated_jobs.append(self)

        def path(self, filename):
            return os.path.join(self.directory, filename)

        def finish(self):
            self.real_time = round(time.time() - self.start_time)
            report.add_job(self)
            name = self.name.lower().replace(" ", "_")
            if name not in report["real_time"]:
                report["real_time"][name] = 0
            self.report["real_time"][name] += job.real_time
            write()

    self.write_report()


    def run(self, args):
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            mlhl = args.unbiased and self.min_rwork > 0.35
            if self.cycle == 1:
                missingPhases = args.colin_hl is None and args.colin_phifom is None
                if missingPhases and args.mr_model is not None:
                    refmac = jobs.Refmac(self.cycle, args.mr_model, mlhl)
                    buccaneer = jobs.Buccaneer(self.cycle, refmac.hklout, args.xyzin, cycles=3)
                else:
                    buccaneer = jobs.Buccaneer(self.cycle, args.mtzin, args.xyzin, cycles=3)
            else:
                coot = jobs.Prune(self.cycle, refmac.xyzout, refmac.hklout)
                refmac = jobs.Refmac(self.cycle, coot.xyzout, mlhl)
                buccaneer = jobs.Buccaneer(self.cycle, refmac.hklout, refmac.xyzout)
            refmac = jobs.Refmac(self.cycle, buccaneer.xyzout, mlhl)
            coot = jobs.Prune(self.cycle, refmac.xyzout, refmac.hklout, chains_only=True)
            if args.add_waters and self.min_rwork < 0.4:
                coot = jobs.FindWaters(self.cycle, coot.xyzout, refmac.hklout)
            refmac = jobs.Refmac(self.cycle, coot.xyzout, mlhl)
            self.process_cycle_output(refmac)
            if self.cycle > 1:
                jobs.remove_job_directories(self.cycle - 1)
            if args.auto_stop and self.cycles_without_improvement == 4:
                break
        jobs.remove_job_directories(self.cycle)

    def improved(self, cycle):
        required_improvement = 0.02
        improvement = (self.min_rwork - cycle["rwork"]) / self.min_rwork
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_built"] - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_sequenced"] - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement:
            return True
        improvement = (self.min_fragments_built - cycle["fragments_built"]) / float(self.min_fragments_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["longest_fragment"] - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement:
            return True
        return False

    def process_cycle_output(self, refmac):
        cycle = {
            "cycle": self.cycle,
            "rwork": refmac.final_rwork,
            "rfree": refmac.final_rfree,
        }
        cycle.update(gemmineer.model_stats(refmac.xyzout))
        print("\nResidues built: %d" % cycle["residues_built"])
        print("Residues sequenced: %d" % cycle["residues_sequenced"])
        print("R-work: %.3f" % cycle["rwork"])
        print("R-free: %.3f" % cycle["rfree"])

        if self.improved(cycle):
            self.cycles_without_improvement = 0
        else:
            self.cycles_without_improvement += 1
            print("\nNo significant improvement for %d cycle(s)" % self.cycles_without_improvement)

        if cycle["rwork"] < self.min_rwork:
            self.min_rwork = cycle["rwork"]
        if cycle["rfree"] < self.min_rfree:
            self.min_rfree = cycle["rfree"]
            print("Copying files to output as R-free improved")
            shutil.copyfile(str(refmac.xyzout), "xyzout.pdb")
            shutil.copyfile(str(refmac.hklout), "hklout.mtz")
        if cycle["residues_built"] > self.max_residues_built:
            self.max_residues_built = cycle["residues_built"]
        if cycle["residues_sequenced"] > self.max_residues_sequenced:
            self.max_residues_sequenced = cycle["residues_sequenced"]
        if cycle["fragments_built"] < self.min_fragments_built:
            self.min_fragments_built = cycle["fragments_built"]
        if cycle["longest_fragment"] > self.max_longest_fragment:
            self.max_longest_fragment = cycle["longest_fragment"]

        report.add_cycle(cycle)

    def remove_job_directories(self, cycle):
        if self.args.keep_intermediate_files:
            return
        for job in self.initiated_jobs:
            if job.cycle == cycle:
                shutil.rmtree(job.directory)
