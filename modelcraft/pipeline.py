from modelcraft.arguments import parse
from modelcraft.buccaneer import Buccaneer
from modelcraft.findwaters import FindWaters
from modelcraft.parrot import Parrot
from modelcraft.prune import Prune
from modelcraft.refmac import Refmac
from modelcraft.sidechains import Sidechains
import json
import shutil
import sys
import time


class Pipeline:
    def __init__(self, argument_list):
        print("# ModelCraft")
        print("\nPlease cite [paper to be published]")
        self.args = parse(argument_list)
        self.initialise()
        self.run()

    def initialise(self):
        self.resolution = self.args.hklin.resolution()
        self.cycle = 0
        self.jobs = {0: []}
        self.current_hkl = self.args.hklin
        self.current_xyz = self.args.xyzin
        self.min_rwork = 100
        self.min_rfree = 100
        self.min_fragments_built = 999
        self.max_longest_fragment = 1
        self.max_residues_built = 1
        self.max_residues_sequenced = 1
        self.cycles_without_improvement = 0
        self.start_time = time.time()
        self.report = {
            "real_time": {"total": 0},
            "cycles": {},
        }

    def run(self):
        args = self.args
        if args.colin_hl is None and args.colin_phifom is None and args.mr_model is not None:
            print("\n## Preparations\n")
            self.get_phases_from_mr_model()
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.jobs[self.cycle] = []
            self.run_cycle()
            self.process_cycle_output()
            self.remove_job_directories(self.cycle - 1)
            if args.auto_stop and self.cycles_without_improvement == 4:
                break
        self.finish()

    def run_cycle(self):
        if self.cycle > 1 and self.resolution < 2.3:
            self.prune()
            self.refmac(cycles=5)
        if self.resolution > 2.4:
            self.parrot()
        self.buccaneer()
        self.refmac(cycles=10)
        self.prune(chains_only=True)
        self.refmac(cycles=5)
        if self.args.fix_side_chains and self.min_rwork < 30:
            self.fix_sidechains()
        if self.min_rwork < 40:
            self.findwaters()

    def finish(self):
        for cycle in range(self.cycle + 1):
            self.remove_job_directories(cycle)
        self.write_report()
        print("\n--- Normal termination ---")
        sys.exit()

    def job_directory(self, name):
        directory = "%02d.%02d_%s" % (self.cycle, len(self.jobs[self.cycle]) + 1, name)
        print(directory)
        return directory

    def add_job(self, job):
        self.jobs[self.cycle].append(job)
        job.finish_time = time.time()
        job.real_time = round(job.finish_time - job.start_time)
        if job.name not in self.report["real_time"]:
            self.report["real_time"][job.name] = 0
        self.report["real_time"][job.name] += job.real_time
        self.write_report()

    def get_phases_from_mr_model(self):
        directory = self.job_directory("mr_refinement")
        job = Refmac(self.args, directory, self.args.mr_model, cycles=10)
        self.add_job(job)
        job.hklout.fphi = None
        self.current_hkl = job.hklout

    def buccaneer(self):
        directory = self.job_directory("buccaneer")
        cycles = 3 if self.cycle == 1 else 2
        job = Buccaneer(self.args, directory, self.current_hkl, self.current_xyz, cycles)
        self.add_job(job)
        if not job.xyzout.exists or job.xyzout.residues == 0:
            print("Stopping the pipeline because buccaneer did not build any residues")
            self.report["cycles"][self.cycle] = {
                "r_work": None,
                "r_free": None,
                "residues": 0,
                "residues_sequenced": 0,
                "waters": 0,
            }
            self.finish()
        self.current_xyz = job.xyzout

    def refmac(self, cycles):
        directory = self.job_directory("refmac")
        use_phases = self.args.unbiased and self.min_rwork > 35
        job = Refmac(self.args, directory, self.current_xyz, cycles, use_phases)
        self.add_job(job)
        self.current_hkl = job.hklout
        self.current_xyz = job.xyzout

    def parrot(self):
        directory = self.job_directory("parrot")
        job = Parrot(self.args, directory, self.current_hkl, self.current_xyz)
        self.add_job(job)
        self.current_hkl = job.hklout

    def prune(self, chains_only=False):
        directory = self.job_directory("prune_chains" if chains_only else "prune")
        job = Prune(directory, self.current_xyz, self.current_hkl, chains_only)
        self.add_job(job)
        self.current_xyz = job.xyzout

    def fix_sidechains(self):
        sidechains_dir = self.job_directory("sidechains")
        sidechains_job = Sidechains(sidechains_dir, self.current_xyz, self.current_hkl)
        self.add_job(sidechains_job)
        refmac_dir = self.job_directory("refmac")
        refmac_job = Refmac(self.args, refmac_dir, sidechains_job.xyzout, 5)
        self.add_job(refmac_job)
        if refmac_job.xyzout.rfree < self.current_xyz.rfree:
            self.current_hkl = refmac_job.hklout
            self.current_xyz = refmac_job.xyzout

    def findwaters(self, dummy=False):
        waters_dir = self.job_directory("finddummys" if dummy else "findwaters")
        waters_job = FindWaters(waters_dir, self.current_xyz, self.current_hkl, dummy)
        self.add_job(waters_job)
        refmac_dir = self.job_directory("refmac")
        refmac_job = Refmac(self.args, refmac_dir, waters_job.xyzout, 10)
        self.add_job(refmac_job)
        if refmac_job.xyzout.rfree < self.current_xyz.rfree:
            self.current_hkl = refmac_job.hklout
            self.current_xyz = refmac_job.xyzout

    def improved(self):
        xyz = self.current_xyz
        required_improvement = 0.02
        improvement = (self.min_rwork - xyz.rwork) / self.min_rwork
        if improvement > required_improvement:
            return True
        improvement = (xyz.residues - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement:
            return True
        improvement = (xyz.sequenced_residues - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement:
            return True
        improvement = (self.min_fragments_built - xyz.fragments) / float(self.min_fragments_built)
        if improvement > required_improvement:
            return True
        improvement = (xyz.longest_fragment - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement:
            return True
        return False

    def process_cycle_output(self):
        print("")
        print("R-work: %.3f" % self.current_xyz.rwork)
        print("R-free: %.3f" % self.current_xyz.rfree)
        print("Residues: %d" % self.current_xyz.residues)
        print("Sequenced residues: %d" % self.current_xyz.sequenced_residues)
        print("Waters: %d" % self.current_xyz.waters)
        self.add_cycle_stats()

        if self.improved():
            self.cycles_without_improvement = 0
        else:
            self.cycles_without_improvement += 1
            print("\nNo significant improvement for %d cycle(s)" % self.cycles_without_improvement)

        if self.current_xyz.rfree < self.min_rfree:
            self.min_rfree = self.current_xyz.rfree
            print("Copying files to output because R-free has improved")
            shutil.copyfile(str(self.current_xyz.path), "modelcraft.pdb")
            shutil.copyfile(str(self.current_hkl.path), "modelcraft.mtz")
            self.add_final_stats()
        self.min_rwork = min(self.min_rwork, self.current_xyz.rwork)
        self.max_residues_built = max(self.max_residues_built, self.current_xyz.residues)
        self.max_residues_sequenced = max(self.max_residues_sequenced, self.current_xyz.sequenced_residues)
        self.min_fragments_built = min(self.min_fragments_built, self.current_xyz.fragments)
        self.max_longest_fragment = max(self.max_longest_fragment, self.current_xyz.longest_fragment)

    def add_cycle_stats(self):
        self.report["cycles"][self.cycle] = self.current_stats()
        self.write_report()

    def add_final_stats(self):
        self.report["final"] = self.current_stats()
        self.write_report()

    def current_stats(self):
        return {
            "r_work": self.current_xyz.rwork,
            "r_free": self.current_xyz.rfree,
            "residues": self.current_xyz.residues,
            "residues_sequenced": self.current_xyz.sequenced_residues,
            "waters": self.current_xyz.waters,
        }

    def write_report(self):
        self.report["real_time"]["total"] = round(time.time() - self.start_time)
        with open("modelcraft.json", "w") as f:
            json.dump(self.report, f, indent=4)

    def remove_job_directories(self, cycle):
        if self.args.keep_intermediate_files:
            return
        for job in self.jobs[cycle]:
            shutil.rmtree(job.directory, ignore_errors=True)
