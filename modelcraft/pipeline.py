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
        self.resolution = self.args.fsigf.resolution_high()
        self.cycle = 0
        self.current_structure = self.args.model
        self.current_phases = self.args.phases
        self.current_fphi_best = None
        self.current_fphi_diff = None
        self.current_fphi_calc = None
        self.last_refmac = None
        self.best_refmac = None
        self.cycles_without_improvement = 0
        self.start_time = time.time()
        self.report = {
            "real_time": {"total": 0},
            "cycles": {},
        }
        self.run()

    def run(self):
        args = self.args
        if args.phases is None and args.model is not None:
            print("\n## Refining Input Model\n")
            self.sheetbend()
            args.model = self.current_structure
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.run_cycle()
            self.process_cycle_output()
            if (
                args.auto_stop
                and self.cycles_without_improvement == args.convergence_cycles
            ):
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
        if self.current_structure is not None:
            self.findwaters(dummy=True)
        self.buccaneer()
        self.prune(chains_only=True)
        self.findwaters()

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---")
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()

    def add_job(self, job):
        if not self.args.keep_jobs:
            job.remove()
        if job.name not in self.report["real_time"]:
            self.report["real_time"][job.name] = 0
        self.report["real_time"][job.name] += job.finish_time - job.start_time
        self.write_report()

    def sheetbend(self):
        print("Sheetbend")
        job = Sheetbend(self.args.fsigf, self.args.freer, self.current_structure)
        self.add_job(job)
        self.refmac(job.structure, cycles=10, auto_accept=True)

    def buccaneer(self):
        print("Buccaneer")
        job = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            input_structure=self.current_structure,
            mr_structure=self.args.model,
            use_mr=True,
            filter_mr=True,
            seed_mr=True,
            cycles=3 if self.cycle == 1 else 2,
            semet=self.args.semet,
            remove_non_protein=self.args.remove_non_protein,
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
        if auto_accept or job.rfree < self.last_refmac.rfree:
            self.update_current_from_refmac_job(job)

    def update_current_from_refmac_job(self, job):
        self.current_structure = job.structure
        self.current_phases = job.abcd
        self.current_fphi_best = job.fphi_best
        self.current_fphi_diff = job.fphi_diff
        self.current_fphi_calc = job.fphi_calc
        self.last_refmac = job

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
        model_stats = ModelStats(self.last_refmac.structure)
        stats = {
            "r_work": self.last_refmac.rwork,
            "r_free": self.last_refmac.rfree,
            "residues": model_stats.residues,
            "sequenced_residues": model_stats.sequenced_residues,
            "fragments": model_stats.fragments,
            "longest_fragment": model_stats.longest_fragment,
            "waters": model_stats.waters,
            "dummy_atoms": model_stats.dummy_atoms,
        }
        self.report["cycles"][self.cycle] = stats
        if self.best_refmac is not None:
            diff = self.best_refmac.rfree - self.last_refmac.rfree
            if diff >= self.args.convergence_tolerance:
                self.cycles_without_improvement = 0
            else:
                self.cycles_without_improvement += 1
        if self.best_refmac is None or self.last_refmac.rfree < self.best_refmac.rfree:
            self.best_refmac = self.last_refmac
            write_mmcif("modelcraft.cif", self.last_refmac.structure)
            write_mtz(
                "modelcraft.mtz",
                [
                    self.args.fsigf,
                    self.args.freer,
                    self.last_refmac.abcd,
                    self.last_refmac.fphi_best,
                    self.last_refmac.fphi_diff,
                    self.last_refmac.fphi_calc,
                ],
            )
            self.report["final"] = stats
        self.write_report()

    def write_report(self):
        self.report["real_time"]["total"] = time.time() - self.start_time
        with open("modelcraft.json", "w") as report_file:
            json.dump(self.report, report_file, indent=4)
