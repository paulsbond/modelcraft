import json
import os
import re
import sys
import time
import gemmi
from . import __version__
from .arguments import parse
from .jobs.buccaneer import Buccaneer
from .jobs.coot import FixSideChains, Prune
from .jobs.ctruncate import CTruncate
from .jobs.findwaters import FindWaters
from .jobs.nautilus import Nautilus
from .jobs.parrot import Parrot
from .jobs.refmac import Refmac, RefmacResult
from .jobs.sheetbend import Sheetbend
from .pipeline import Pipeline
from .reflections import DataItem, write_mtz
from .structure import ModelStats, write_mmcif


class ModelCraft(Pipeline):
    def __init__(self, argument_list):
        print(f"# ModelCraft {__version__}\n")
        print("Arguments:")
        print(" %s\n" % " ".join(argument_list).replace(" --", "\n --"))
        self.args = parse(argument_list)
        super().__init__(keep_jobs=self.args.keep_jobs, keep_logs=self.args.keep_logs)
        self.cycle = 0
        self.current_structure: gemmi.Structure = self.args.model
        self.current_phases: DataItem = self.args.phases
        self.current_fphi_best: DataItem = None
        self.current_fphi_diff: DataItem = None
        self.current_fphi_calc: DataItem = None
        self.last_refmac: RefmacResult = None
        self.best_refmac: RefmacResult = None
        self.cycles_without_improvement = 0
        self.start_time = None
        self.report = {
            "seconds": self.seconds,
            "cycles": {},
        }

    @property
    def resolution(self):
        return self.args.fsigf.resolution_high()

    def run(self):
        self.start_time = time.time()
        args = self.args
        os.makedirs(args.directory, exist_ok=True)
        os.chdir(args.directory)
        _check_for_files_that_could_be_overwritten()
        if self.args.observations.types == "FQ":
            self.args.fsigf = self.args.observations
        else:
            print("\n## Converting input observations to mean amplitudes\n")
            result = CTruncate(observations=self.args.observations).run(self)
            self.args.fsigf = result.fmean
        if args.model is not None:
            print("\n## Refining Input Model\n")
            self.sheetbend()
            args.model = self.current_structure
            if args.phases is not None:
                self.current_phases = args.phases
            _print_refmac_result(self.last_refmac)
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.run_cycle()
            self.process_cycle_output()
            if (
                args.auto_stop
                and self.cycles_without_improvement == args.convergence_cycles
            ):
                break
        if not args.basic and self.best_refmac.rwork < 30 and self.resolution < 2.5:
            print("\n## Finalisations\n")
            self.cycle += 1
            self.update_current_from_refmac_result(self.best_refmac)
            self.fixsidechains()
            self.process_cycle_output()
        print("\n## Best Model:")
        _print_refmac_result(self.best_refmac)
        self.terminate(reason="Normal")

    def run_cycle(self):
        if self.args.basic:
            if self.cycle == 1:
                self.parrot()
            self.buccaneer()
            self.nautilus()
        else:
            if self.cycle > 1 and self.resolution < 2.3:
                self.prune()
            self.parrot()
            if self.current_structure is not None:
                self.findwaters(dummy=True)
            self.buccaneer()
            self.prune(chains_only=True)
            self.nautilus()
            self.findwaters()

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---")
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()

    def sheetbend(self):
        print("Sheetbend")
        result = Sheetbend(
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            structure=self.current_structure,
            executable=self.args.sheetbend,
        ).run(self)
        self.refmac(result.structure, cycles=10, auto_accept=True)

    def buccaneer(self):
        if not self.args.contents.proteins:
            return
        print("Buccaneer")
        result = Buccaneer(
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
            executable=self.args.buccaneer,
        ).run(self)
        stats = ModelStats(result.structure)
        if stats.residues == 0:
            self.terminate(reason="Buccaneer did not build any residues")
        self.refmac(result.structure, cycles=10, auto_accept=True)

    def nautilus(self):
        if not (self.args.contents.rnas or self.args.contents.dnas):
            return
        print("Nautilus")
        result = Nautilus(
            contents=self.args.contents,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            structure=self.current_structure,
        ).run(self)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def refmac(self, structure: gemmi.Structure, cycles: int, auto_accept: bool):
        print("REFMAC")
        use_phases = self.args.unbiased and (
            self.best_refmac is None or self.best_refmac.rwork > 35
        )
        result = Refmac(
            structure=structure,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            cycles=cycles,
            phases=self.args.phases if use_phases else None,
            twinned=self.args.twinned,
        ).run(self)
        if auto_accept or result.rfree < self.last_refmac.rfree:
            self.update_current_from_refmac_result(result)

    def update_current_from_refmac_result(self, result: RefmacResult):
        self.current_structure = result.structure
        self.current_phases = result.abcd
        self.current_fphi_best = result.fphi_best
        self.current_fphi_diff = result.fphi_diff
        self.current_fphi_calc = result.fphi_calc
        self.last_refmac = result

    def parrot(self):
        print("Parrot")
        result = Parrot(
            contents=self.args.contents,
            fsigf=self.args.fsigf,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            structure=self.current_structure,
            executable=self.args.parrot,
        ).run(self)
        self.current_phases = result.abcd
        self.current_fphi_best = result.fphi

    def prune(self, chains_only=False):
        if not self.args.contents.proteins:
            return
        print("Pruning chains" if chains_only else "Pruning model")
        result = Prune(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
            chains_only=chains_only,
        ).run(self)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def fixsidechains(self):
        if not self.args.contents.proteins:
            return
        print("Fixing side chains")
        result = FixSideChains(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
        ).run(self)
        self.refmac(result.structure, cycles=5, auto_accept=False)

    def findwaters(self, dummy=False):
        print("Adding dummy atoms" if dummy else "Adding waters")
        result = FindWaters(
            structure=self.current_structure,
            fphi=self.current_fphi_best,
            dummy=dummy,
        ).run(self)
        self.refmac(result.structure, cycles=10, auto_accept=False)

    def process_cycle_output(self):
        _print_refmac_result(self.last_refmac)
        model_stats = ModelStats(self.last_refmac.structure)
        stats = {
            "residues": model_stats.residues,
            "waters": model_stats.waters,
            "r_work": self.last_refmac.rwork,
            "r_free": self.last_refmac.rfree,
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
        self.seconds["total"] = time.time() - self.start_time
        with open("modelcraft.json", "w") as report_file:
            json.dump(self.report, report_file, indent=4)


def _print_refmac_result(result: RefmacResult):
    model_stats = ModelStats(result.structure)
    print("")
    print(f"Residues: {model_stats.residues:5d}")
    print(f"Waters:   {model_stats.waters:5d}")
    print(f"R-work:   {result.rwork:5.1f}")
    print(f"R-free:   {result.rfree:5.1f}")


def _check_for_files_that_could_be_overwritten():
    patterns = [
        r"modelcraft\.cif",
        r"modelcraft\.json",
        r"modelcraft\.mtz",
        r"job_[A-Za-z0-9]+_[A-Za-z0-9]{20}",
        r"job_\d+_[A-Za-z0-9]+",
    ]
    paths = os.listdir(".")
    paths = [p for p in paths if any(re.fullmatch(pattern, p) for pattern in patterns)]
    if paths:
        print("\nThe following files may be from a previous run:\n")
        for path in paths:
            print("-", path)
        print("\nPlease run in a different directory or remove these files.")
        sys.exit()
