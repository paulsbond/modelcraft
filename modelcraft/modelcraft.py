import dataclasses
import json
import os
import sys
import time
import gemmi
from . import __version__
from .arguments import parse
from .jobs.buccaneer import Buccaneer, BuccaneerResult
from .jobs.coot import FixSideChains, Prune
from .jobs.ctruncate import CTruncate
from .jobs.findwaters import FindWaters, FindWatersResult
from .jobs.nautilus import Nautilus, NautilusResult
from .jobs.parrot import Parrot, ParrotResult
from .jobs.refmac import RefmacXray, RefmacEm, RefmacResult
from .jobs.sheetbend import Sheetbend, SheetbendResult
from .cell import max_distortion, remove_scale, update_cell
from .pipeline import Pipeline
from .reflections import DataItem, write_mtz
from .structure import ModelStats, remove_residues, write_mmcif


class ModelCraft(Pipeline):
    def __init__(self, args):
        self.args = parse(args)
        super().__init__(
            directory=self.args.directory,
            keep_jobs=self.args.keep_files,
            keep_logs=self.args.keep_logs,
        )
        self.cycle = 0

        self.output_refmac: RefmacResult = None
        self.cycles_without_improvement = 0
        self.start_time = None
        self.report = {
            "version": __version__,
            "args": args or sys.argv[1:],
            "seconds": self.seconds,
            "cycles": [],
            "jobs": [],
        }

    @property
    def resolution(self):
        return self.args.fmean.resolution_high()

    def run(self):
        self.start_time = time.time()
        args = self.args
        print(f"# ModelCraft {__version__}")
        os.makedirs(args.directory, exist_ok=False)
        self.convert_observations()
        # sheetbend = self.sheetbend(self.args.model)
        refmac = self.refmac(self.args.model, cycles=5)
        parrot = self.parrot(refmac.abcd, refmac.fphi_best)
        findwaters = self.findwaters(refmac.structure, parrot.fphi, dummy=True)
        refmac = self.refmac(findwaters.structure, cycles=5)
        remove_residues(refmac.structure, "DUM")
        nautilus = self.nautilus(refmac.abcd, refmac.fphi_best, refmac.structure)
        # buccaneer = self.buccaneer(refmac.abcd, refmac.fphi_best, refmac.structure)
        self.terminate(reason="Normal")

    def convert_observations(self):
        if self.args.fmean is None:
            observations = self.args.ianom or self.args.imean or self.args.fanom
            self._running_job("CTruncate")
            ctruncate = CTruncate(observations=observations).run(self)
            self._finished_job("CTruncate", ctruncate)
            self.args.fmean = ctruncate.fmean
            if self.args.fanom is None and ctruncate.fanom is not None:
                self.args.fanom = ctruncate.fanom
            if self.args.imean is None and ctruncate.imean is not None:
                self.args.imean = ctruncate.imean

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---")
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()

    def sheetbend(self, structure: gemmi.Structure) -> SheetbendResult:
        self._running_job("Sheetbend")
        result = Sheetbend(structure, self.args.fmean, self.args.freer).run(self)
        self._finished_job("Sheetbend", result)
        return result

    def buccaneer(
        self, phases: DataItem, fphi: DataItem = None, structure: gemmi.Structure = None
    ) -> BuccaneerResult:
        if not self.args.contents.proteins:
            return
        self._running_job("Buccaneer")
        result = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=phases,
            fphi=fphi,
            freer=self.args.freer,
            input_structure=structure,
            mr_structure=self.args.model,
            cycles=2,
            em_mode=self.args.mode == "em",
        ).run(self)
        self._finished_job("Buccaneer", result)
        return result

    def nautilus(
        self, phases: DataItem, fphi: DataItem = None, structure: gemmi.Structure = None
    ) -> NautilusResult:
        if not (self.args.contents.rnas or self.args.contents.dnas):
            return
        self._running_job("Nautilus")
        result = Nautilus(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=phases,
            fphi=fphi,
            freer=self.args.freer,
            structure=structure,
            cycles=1,
        ).run(self)
        self._finished_job("Nautilus", result)
        return result

    def refmac(self, structure: gemmi.Structure, cycles: int) -> RefmacResult:
        self._running_job("Refmac")
        if self.args.mode == "xray":
            use_phases = self.args.unbiased and (
                self.output_refmac is None or self.output_refmac.rwork > 0.35
            )
            result = RefmacXray(
                structure=structure,
                fsigf=self.args.fmean,
                freer=self.args.freer,
                cycles=cycles,
                phases=self.args.phases if use_phases else None,
                twinned=self.args.twinned,
            ).run(self)
        else:
            result = RefmacEm(
                structure=structure,
                fphi=self.args.fphi,
                cycles=cycles,
            ).run(self)
        self._finished_job("Refmac", result)
        return result

    def parrot(
        self, phases: DataItem, fphi: DataItem = None, structure: gemmi.Structure = None
    ) -> ParrotResult:
        self._running_job("Parrot")
        result = Parrot(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            freer=self.args.freer,
            phases=phases,
            fphi=fphi,
            structure=structure,
        ).run(self)
        self._finished_job("Parrot", result)
        return result

    def prune(self, chains_only=False):
        if self.args.disable_pruning or not self.args.contents.proteins:
            return
        name = "Pruning chains" if chains_only else "Pruning model"
        self._running_job(name)
        result = Prune(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
            chains_only=chains_only,
        ).run(self)
        self._finished_job(name, result)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def fixsidechains(self):
        if self.args.disable_side_chain_fixing or not self.args.contents.proteins:
            return
        self._running_job("Fixing side chains")
        result = FixSideChains(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
        ).run(self)
        self._finished_job("Fixing side chains", result)
        self.refmac(result.structure, cycles=5, auto_accept=False)

    def findwaters(
        self, structure: gemmi.Structure, fphi: DataItem, dummy=False
    ) -> FindWatersResult:
        name = "Adding dummy atoms" if dummy else "Adding waters"
        self._running_job(name)
        result = FindWaters(structure=structure, fphi=fphi, dummy=dummy).run(self)
        self._finished_job(name, result)
        return result

    def process_cycle_output(self, result: RefmacResult):
        self.print_refmac_result(result)
        model_stats = ModelStats(result.structure)
        stats = {"cycle": self.cycle, "residues": model_stats.residues}
        if self.args.mode == "xray":
            stats["waters"] = model_stats.waters
            stats["r_work"] = result.rwork
            stats["r_free"] = result.rfree
        if self.args.mode == "em":
            stats["fsc"] = result.fsc
        self.report["cycles"].append(stats)
        if (
            self.output_refmac is None
            or (self.args.mode == "xray" and result.rwork < self.output_refmac.rwork)
            or (self.args.mode == "em" and result.fsc > self.output_refmac.fsc)
        ):
            self.cycles_without_improvement = 0
            self.output_refmac = result
            write_mmcif(self.path("modelcraft.cif"), result.structure)
            write_mtz(
                self.path("modelcraft.mtz"),
                [
                    self.args.fmean,
                    self.args.freer,
                    result.abcd,
                    result.fphi_best,
                    result.fphi_diff,
                    result.fphi_calc,
                ],
            )
            self.report["final"] = stats
        else:
            self.cycles_without_improvement += 1
        self.write_report()

    def _running_job(self, name):
        print(name)
        self.report["running_job"] = name
        self.write_report()

    def _finished_job(self, name, result):
        self.report.pop("running_job", None)
        result_dict = {}
        for field in dataclasses.fields(result):
            value = getattr(result, field.name)
            try:
                json.dumps(value)
            except TypeError:
                pass
            else:
                result_dict[field.name] = value
        print(json.dumps(result_dict, indent=4))
        self.report["jobs"].append({"name": name, **result_dict})
        self.write_report()

    def write_report(self):
        self.seconds["total"] = time.time() - self.start_time
        with open(self.path("modelcraft.json"), "w") as report_file:
            json.dump(self.report, report_file, indent=4)

    def print_refmac_result(self, result: RefmacResult):
        model_stats = ModelStats(result.structure)
        print(f"\nResidues: {model_stats.residues:6d}")
        if self.args.mode == "xray":
            print(f"Waters:   {model_stats.waters:6d}")
            print(f"R-work:   {result.rwork:6.4f}")
            print(f"R-free:   {result.rfree:6.4f}")
        if self.args.mode == "em":
            print(f"FSC:      {result.fsc:6.4f}")

    def update_model_cell(self):
        structure = self.args.model
        mtz = self.args.fmean
        structure_spacegroup = gemmi.find_spacegroup_by_name(
            structure.spacegroup_hm,
            alpha=structure.cell.alpha,
            gamma=structure.cell.gamma,
        )
        distortion = max_distortion(old_cell=structure.cell, new_cell=mtz.cell)
        if structure_spacegroup.number != mtz.spacegroup.number or distortion > 0.05:
            print("The model cell is incompatible with the data cell")
            cell1 = " ".join(f"{x:7.2f}" for x in structure.cell.parameters)
            cell2 = " ".join(f"{x:7.2f}" for x in mtz.cell.parameters)
            print(f"Model: {cell1}  {structure_spacegroup.hm}")
            print(f"Data:  {cell2}  {mtz.spacegroup.hm}")
            print("Molecular replacement should be used first")
            self.terminate("Model cell is incompatible")
        remove_scale(structure=structure)
        if distortion > 0:
            update_cell(structure=structure, new_cell=mtz.cell)
