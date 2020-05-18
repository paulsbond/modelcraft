import distutils.spawn
import json
import sys
import time
import gemmi
from .arguments import parse
from .jobs import (
    Buccaneer,
    FindWaters,
    FixSideChains,
    Parrot,
    Prune,
    Refmac,
    Sheetbend,
)
from .reflections import write_mtz
from .structure import ModelStats, write_mmcif


class Pipeline:
    def __init__(self, argument_list):
        print("# ModelCraft\n")
        print("Arguments:")
        print(" %s\n" % " ".join(argument_list).replace(" --", "\n --"))
        self.args = parse(argument_list)
        self.resolution = self.args.fsigf.resolution  # TODO: Change after gemmi update
        # self.resolution = self.args.fsigf.resolution_high()  # To this
        self.cycle = 0
        self.current_structure = self.args.xyzin
        self.current_phases = self.args.phases
        self.current_fphi_best = None
        self.current_fphi_diff = None
        self.current_fphi_calc = None
        self.current_rwork = None
        self.current_rfree = None
        self.best_refmac = None
        self.best_refmac_cycle = 0
        self.start_time = time.time()
        self.report = {
            "real_time": {"total": 0},
            "cycles": {},
            "termination_reason": "Unknown",
        }
        self.run()

    def run(self):
        args = self.args
        if args.phases is None and args.mr_model is not None:
            print("\n## Preparations\n")
            self.get_phases_from_mr_model()
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.run_cycle()
            self.process_cycle_output()
            if args.auto_stop and self.cycle - self.best_refmac_cycle == 4:
                break
        if self.best_refmac.rwork < 30 and self.resolution < 2.5:
            print("\n## Finalisations\n")
            self.cycle += 1
            self.update_current_from_refmac_job(self.best_refmac)
            self.fixsidechains()
            self.process_cycle_output()
        self.terminate(reason="Normal")

    def run_cycle(self):
        if self.cycle > 1 and self.resolution < 2.3:
            self.prune()
        self.parrot()
        if self.current_structure is not None and self.current_fphi_best is not None:
            self.findwaters(dummy=True)
        self.buccaneer()
        self.prune(chains_only=True)
        if self.best_refmac.rwork < 40:
            self.findwaters()

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---")
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()

    def add_job(self, job):
        if job.name not in self.report["real_time"]:
            self.report["real_time"][job.name] = 0
        self.report["real_time"][job.name] += job.finish_time - job.start_time
        self.write_report()

    def get_phases_from_mr_model(self):
        structure = self.args.mr_model
        if distutils.spawn.find_executable("csheetbend") is not None:
            print("Sheetbend")
            sheetbend = Sheetbend(self.args.fsigf, self.args.freer, structure)
            self.add_job(sheetbend)
            structure = sheetbend.structure
        self.refmac(structure, cycles=10, auto_accept=True)

    def buccaneer(self):
        print("Buccaneer")
        job = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            input_structure=self.current_structure,
            known_structure=self.args.known_structure,
            mr_structure=self.args.mr_model if self.args.mr_mode > 1 else None,
            use_mr=self.args.mr_mode > 2,
            filter_mr=self.args.mr_mode in (4, 6),
            seed_mr=self.args.mr_mode > 4,
            cycles=3 if self.cycle == 1 else 2,
            semet=self.args.semet,
            program=self.args.buccaneer,
        )
        self.add_job(job)
        stats = ModelStats(job.structure)
        if stats.residues == 0:
            self.terminate(reason="Buccaneer did not build any residues")
        self.refmac(job.structure, cycles=10, auto_accept=True)

    def refmac(self, structure: gemmi.Structure, cycles: int, auto_accept: bool):
        print("REFMAC")
        use_phases = self.args.unbiased and (
            self.best_refmac is None or self.best_refmac.rwork > 35
        )
        job = Refmac(
            structure=structure,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            cycles=cycles,
            phases=self.args.phases if use_phases else None,
            twinned=self.args.twinned,
        )
        self.add_job(job)
        if auto_accept or job.rfree < self.current_rfree:
            self.update_current_from_refmac_job(job)
            if self.best_refmac is None or job.rfree < self.best_refmac.rfree:
                self.best_refmac = job
                self.best_refmac_cycle = self.cycle
                write_mmcif("modelcraft.cif", job.structure)
                write_mtz(
                    "modelcraft.mtz",
                    [
                        self.args.fsigf,
                        self.args.freer,
                        job.abcd,
                        job.fphi_best,
                        job.fphi_diff,
                        job.fphi_calc,
                    ],
                )
                self.report["final"] = self.current_stats()
                self.write_report()

    def update_current_from_refmac_job(self, job):
        self.current_structure = job.structure
        self.current_phases = job.abcd
        self.current_fphi_best = job.fphi_best
        self.current_fphi_diff = job.fphi_diff
        self.current_fphi_calc = job.fphi_calc
        self.current_rwork = job.rwork
        self.current_rfree = job.rfree

    def parrot(self):
        print("Parrot")
        job = Parrot(
            contents=self.args.contents,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            structure=self.current_structure,
        )
        self.add_job(job)
        self.current_phases = job.abcd
        self.current_fphi_best = job.fphi

    def prune(self, chains_only=False):
        print("Pruning chains" if chains_only else "Pruning model")
        job = Prune(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
            chains_only=chains_only,
        )
        self.add_job(job)
        self.refmac(job.structure, cycles=5, auto_accept=True)

    def fixsidechains(self):
        print("Fixing side chains")
        job = FixSideChains(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
        )
        self.add_job(job)
        self.refmac(job.structure, cycles=5, auto_accept=False)

    def findwaters(self, dummy=False):
        print("Adding dummy atoms" if dummy else "Adding waters")
        job = FindWaters(
            structure=self.current_structure, fphi=self.current_fphi_best, dummy=dummy
        )
        self.add_job(job)
        self.refmac(job.structure, cycles=10, auto_accept=False)

    def process_cycle_output(self):
        self.report["cycles"][self.cycle] = self.current_stats()
        self.write_report()

    def current_stats(self):
        stats = ModelStats(self.current_structure)
        return {
            "r_work": self.current_rwork,
            "r_free": self.current_rfree,
            "residues": stats.residues,
            "sequenced_residues": stats.sequenced_residues,
            "fragments": stats.fragments,
            "longest_fragment": stats.longest_fragment,
            "waters": stats.waters,
            "dummy_atoms": stats.dummy_atoms,
        }

    def write_report(self):
        self.report["real_time"]["total"] = time.time() - self.start_time
        with open("modelcraft.json", "w") as report_file:
            json.dump(self.report, report_file, indent=4)
