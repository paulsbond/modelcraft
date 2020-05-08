import distutils.spawn
import json
import os
import shutil
import sys
import time
from modelcraft.arguments import parse
from modelcraft.jobs import (
    Buccaneer,
    FindWaters,
    FixSideChains,
    Parrot,
    Prune,
    Refmac,
    Sheetbend,
)


class Pipeline:
    def __init__(self, argument_list):
        print("# ModelCraft\n")
        print("Arguments:")
        print(" %s\n" % " ".join(argument_list).replace(" --", "\n --"))
        self.args = parse(argument_list)
        self.resolution = self.args.fsigf
        self.cycle = 0
        self.current_hkl = self.args.hklin
        self.current_xyz = self.args.xyzin
        self.best_refmac: Refmac = None
        self.cycles_without_improvement = 0
        self.start_time = time.time()
        self.report = {
            "real_time": {"total": 0},
            "cycles": {},
        }
        self.run()

    def run(self):
        args = self.args
        if args.phases is None and args.mr_model is not None:
            print("\n## Preparations\n")
    #         self.get_phases_from_mr_model()
    #     for self.cycle in range(1, args.cycles + 1):
    #         print("\n## Cycle %d\n" % self.cycle)
    #         self.run_cycle()
    #         self.process_cycle_output()
    #         if args.auto_stop and self.cycles_without_improvement == 4:
    #             break
    #     if self.min_rwork < 30 and self.resolution < 2.5:
    #         print("\n## Finalisations\n")
    #         self.cycle += 1
    #         self.current_xyz = self.best_xyz
    #         self.current_hkl = self.best_hkl
    #         self.fixsidechains()
    #         self.process_cycle_output()
    #     self.finish()

    # def run_cycle(self):
    #     if self.cycle > 1 and self.resolution < 2.3:
    #         self.prune()
    #         self.refmac(cycles=5)
    #     if self.resolution > 2.4:
    #         self.parrot()
    #     self.buccaneer()
    #     self.refmac(cycles=10)
    #     self.prune(chains_only=True)
    #     self.refmac(cycles=5)
    #     if self.min_rwork < 40:
    #         self.findwaters()

    # def finish(self):
    #     self.write_report()
    #     print("\n--- Normal termination ---")
    #     sys.exit()

    # def add_job(self, job):
    #     if job.name not in self.report["real_time"]:
    #         self.report["real_time"][job.name] = 0
    #     self.report["real_time"][job.name] += job.finish_time - job.start_time
    #     self.write_report()

    # def get_phases_from_mr_model(self):
    #     xyzin = self.args.mr_model
    #     if distutils.spawn.find_executable("csheetbend") is not None:
    #         sheetbend_job = Sheetbend(self.args, self.args.mr_model)
    #         self.add_job(sheetbend_job)
    #         xyzin = sheetbend_job.xyzout
    #     refmac_job = Refmac(self.args, xyzin, cycles=10)
    #     self.add_job(refmac_job)
    #     self.current_hkl = refmac_job.hklout

    # def buccaneer(self):
    #     cycles = 3 if self.cycle == 1 else 2
    #     job = Buccaneer(self.args, self.current_hkl, self.current_xyz, cycles)
    #     self.add_job(job)
    #     if not job.xyzout.exists or job.xyzout.residues == 0:
    #         print("Stopping the pipeline because buccaneer did not build any residues")
    #         self.report["cycles"][self.cycle] = {
    #             "r_work": None,
    #             "r_free": None,
    #             "residues": 0,
    #             "residues_sequenced": 0,
    #             "waters": 0,
    #         }
    #         self.finish()
    #     self.current_xyz = job.xyzout

    # def refmac(self, cycles):
    #     use_phases = self.args.unbiased and self.min_rwork > 35
    #     job = Refmac(self.args, self.current_xyz, cycles, use_phases)
    #     self.add_job(job)
    #     self.current_hkl = job.hklout
    #     self.current_xyz = job.xyzout

    # def parrot(self):
    #     job = Parrot(self.args, self.current_hkl, self.current_xyz)
    #     self.add_job(job)
    #     self.current_hkl = job.hklout

    # def prune(self, chains_only=False):
    #     job = Prune(self.current_xyz, self.current_hkl, chains_only)
    #     self.add_job(job)
    #     self.current_xyz = job.xyzout

    # def fixsidechains(self):
    #     sidechains_job = FixSidechains(self.current_xyz, self.current_hkl)
    #     self.add_job(sidechains_job)
    #     refmac_job = Refmac(self.args, refmac_dir, sidechains_job.xyzout, 5)
    #     self.add_job(refmac_job)
    #     if refmac_job.xyzout.rfree < self.current_xyz.rfree:
    #         self.current_hkl = refmac_job.hklout
    #         self.current_xyz = refmac_job.xyzout

    # def findwaters(self, dummy=False):
    #     waters_job = FindWaters(self.current_xyz, self.current_hkl, dummy)
    #     self.add_job(waters_job)
    #     refmac_job = Refmac(self.args, waters_job.xyzout, 10)
    #     self.add_job(refmac_job)
    #     if refmac_job.xyzout.rfree < self.current_xyz.rfree:
    #         self.current_hkl = refmac_job.hklout
    #         self.current_xyz = refmac_job.xyzout

    # def improved(self):
    #     xyz = self.current_xyz
    #     required_improvement = 0.02
    #     improvement = (self.min_rwork - xyz.rwork) / self.min_rwork
    #     if improvement > required_improvement:
    #         return True
    #     improvement = (xyz.residues - self.max_residues_built) / float(
    #         self.max_residues_built
    #     )
    #     if improvement > required_improvement:
    #         return True
    #     improvement = (xyz.sequenced_residues - self.max_residues_sequenced) / float(
    #         self.max_residues_sequenced
    #     )
    #     if improvement > required_improvement:
    #         return True
    #     improvement = (self.min_fragments_built - xyz.fragments) / float(
    #         self.min_fragments_built
    #     )
    #     if improvement > required_improvement:
    #         return True
    #     improvement = (xyz.longest_fragment - self.max_longest_fragment) / float(
    #         self.max_longest_fragment
    #     )
    #     if improvement > required_improvement:
    #         return True
    #     return False

    # def process_cycle_output(self):
    #     print("")
    #     print("R-work: %.3f" % self.current_xyz.rwork)
    #     print("R-free: %.3f" % self.current_xyz.rfree)
    #     print("Residues: %d" % self.current_xyz.residues)
    #     print("Sequenced residues: %d" % self.current_xyz.sequenced_residues)
    #     print("Waters: %d" % self.current_xyz.waters)
    #     self.add_cycle_stats()

    #     if self.improved():
    #         self.cycles_without_improvement = 0
    #     else:
    #         self.cycles_without_improvement += 1
    #         print(
    #             "\nNo significant improvement for %d cycle(s)"
    #             % self.cycles_without_improvement
    #         )

    #     if self.current_xyz.rfree < self.min_rfree:
    #         self.min_rfree = self.current_xyz.rfree
    #         print("Copying files to output because R-free has improved")
    #         shutil.copyfile(str(self.current_xyz.path), "modelcraft.pdb")
    #         shutil.copyfile(str(self.current_hkl.path), "modelcraft.mtz")
    #         self.best_xyz = self.current_xyz
    #         self.best_hkl = self.current_hkl
    #         self.best_xyz.path = os.path.abspath("modelcraft.pdb")
    #         self.best_hkl.path = os.path.abspath("modelcraft.mtz")
    #         self.add_final_stats()
    #     self.min_rwork = min(self.min_rwork, self.current_xyz.rwork)
    #     self.max_residues_built = max(
    #         self.max_residues_built, self.current_xyz.residues
    #     )
    #     self.max_residues_sequenced = max(
    #         self.max_residues_sequenced, self.current_xyz.sequenced_residues
    #     )
    #     self.min_fragments_built = min(
    #         self.min_fragments_built, self.current_xyz.fragments
    #     )
    #     self.max_longest_fragment = max(
    #         self.max_longest_fragment, self.current_xyz.longest_fragment
    #     )

    # def add_cycle_stats(self):
    #     self.report["cycles"][self.cycle] = self.current_stats()
    #     self.write_report()

    # def add_final_stats(self):
    #     self.report["final"] = self.current_stats()
    #     self.write_report()

    # def current_stats(self):
    #     return {
    #         "r_work": self.current_xyz.rwork,
    #         "r_free": self.current_xyz.rfree,
    #         "residues": self.current_xyz.residues,
    #         "residues_sequenced": self.current_xyz.sequenced_residues,
    #         "waters": self.current_xyz.waters,
    #     }

    # def write_report(self):
    #     self.report["real_time"]["total"] = round(time.time() - self.start_time)
    #     with open("modelcraft.json", "w") as f:
    #         json.dump(self.report, f, indent=4)
