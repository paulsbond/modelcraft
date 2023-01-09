import dataclasses
import json
import os
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
from .jobs.refmac import RefmacXray, RefmacEm, RefmacResult
from .jobs.sheetbend import Sheetbend
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
        self.current_structure: gemmi.Structure = self.args.model
        self.current_phases: DataItem = self.args.phases
        self.current_fphi_best: DataItem = None
        self.current_fphi_diff: DataItem = None
        self.current_fphi_calc: DataItem = None
        self.last_refmac: RefmacResult = None
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
        print(f"# ModelCraft {__version__}", flush=True)
        os.makedirs(args.directory, exist_ok=False)
        self._convert_observations()
        self._refine_input_model()
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle, flush=True)
            self.run_cycle()
            self.process_cycle_output(self.last_refmac)
            if (
                self.cycles_without_improvement == args.auto_stop_cycles
                and args.auto_stop_cycles > 0
            ):
                break
        if (
            args.mode == "xray"
            and not args.basic
            and self.output_refmac.rwork < 0.3
            and self.resolution < 2.5
        ):
            print("\n## Finalisations\n", flush=True)
            self.cycle += 1
            self.update_current_from_refmac_result(self.output_refmac)
            self.fixsidechains()
            self.process_cycle_output(self.last_refmac)
        print("\n## Best Model:", flush=True)
        self.print_refmac_result(self.output_refmac)
        self.terminate(reason="Normal")

    def _convert_observations(self):
        if self.args.fmean is None:
            print("\n## Converting input observations to mean amplitudes\n", flush=True)
            observations = self.args.ianom or self.args.imean or self.args.fanom
            self._running_job("CTruncate")
            ctruncate = CTruncate(observations=observations).run(self)
            self._finished_job("CTruncate", ctruncate)
            self.args.fmean = ctruncate.fmean
            if self.args.fanom is None and ctruncate.fanom is not None:
                self.args.fanom = ctruncate.fanom
            if self.args.imean is None and ctruncate.imean is not None:
                self.args.imean = ctruncate.imean

    def _refine_input_model(self):
        if self.args.mode == "xray" and self.args.model is not None:
            print("\n## Refining Input Model\n", flush=True)
            self.update_model_cell()
            self.sheetbend()
            self.args.model = self.current_structure
            if self.args.phases is not None:
                self.current_phases = self.args.phases
            self.print_refmac_result(self.last_refmac)

    def run_cycle(self):
        if self.args.mode == "em":
            self.buccaneer()
            self.nautilus()
        elif self.args.basic:
            if self.cycle == 1:
                self.parrot()
            self.buccaneer()
            self.nautilus()
        else:
            if self.cycle > 1 and self.resolution < 2.3:
                self.prune()
            self.parrot()
            if self.current_structure is not None:
                if self.cycle > 1 or self.args.phases is None:
                    self.findwaters(dummy=True)
                remove_residues(structure=self.current_structure, names={"HOH", "DUM"})
            self.buccaneer()
            self.prune(chains_only=True)
            self.nautilus()
            self.findwaters()

    def terminate(self, reason: str):
        print(f"\n--- Termination: {reason} ---", flush=True)
        self.report["termination_reason"] = reason
        self.write_report()
        sys.exit()

    def sheetbend(self):
        if self.args.disable_sheetbend:
            self.refmac(self.current_structure, cycles=10, auto_accept=True)
        else:
            self._running_job("Sheetbend")
            result = Sheetbend(
                fsigf=self.args.fmean,
                freer=self.args.freer,
                structure=self.current_structure,
            ).run(self)
            self._finished_job("Sheetbend", result)
            self.refmac(result.structure, cycles=10, auto_accept=True)

    def buccaneer(self):
        if not self.args.contents.proteins:
            return
        self._running_job("Buccaneer")
        result = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.current_phases,
            fphi=self.current_fphi_best if self.args.mode == "xray" else None,
            freer=self.args.freer if self.args.mode == "xray" else None,
            input_structure=self.current_structure,
            mr_structure=self.args.model,
            use_mr=True,
            filter_mr=True,
            seed_mr=True,
            cycles=3 if self.cycle == 1 else 2,
            em_mode=self.args.mode == "em",
        ).run(self)
        self._finished_job("Buccaneer", result)
        if result.residues_built == 0:
            self.terminate(reason="Buccaneer did not build any residues")
        self.refmac(result.structure, cycles=10, auto_accept=True)

    def nautilus(self):
        if not (self.args.contents.rnas or self.args.contents.dnas):
            return
        self._running_job("Nautilus")
        result = Nautilus(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.current_phases,
            fphi=self.current_fphi_best if self.args.mode == "xray" else None,
            freer=self.args.freer if self.args.mode == "xray" else None,
            structure=self.current_structure,
        ).run(self)
        self._finished_job("Nautilus", result)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def refmac(self, structure: gemmi.Structure, cycles: int, auto_accept: bool):
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
        if auto_accept or result.rfree < self.last_refmac.rfree:
            if not auto_accept:
                print("(accepted)", flush=True)
            self.update_current_from_refmac_result(result)
        else:
            print("(rejected)", flush=True)

    def update_current_from_refmac_result(self, result: RefmacResult):
        self.current_structure = result.structure
        self.current_phases = result.abcd
        self.current_fphi_best = result.fphi_best
        self.current_fphi_diff = result.fphi_diff
        self.current_fphi_calc = result.fphi_calc
        self.last_refmac = result

    def parrot(self):
        if self.args.disable_parrot:
            return
        self._running_job("Parrot")
        result = Parrot(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            structure=self.current_structure,
        ).run(self)
        self._finished_job("Parrot", result)
        self.current_phases = result.abcd
        self.current_fphi_best = result.fphi

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

    def findwaters(self, dummy=False):
        if dummy and self.args.disable_dummy_atoms:
            return
        if not dummy and self.args.disable_waters:
            return
        name = "Adding dummy atoms" if dummy else "Adding waters"
        self._running_job(name)
        result = FindWaters(
            structure=self.current_structure,
            fphi=self.current_fphi_best,
            dummy=dummy,
        ).run(self)
        self._finished_job(name, result)
        self.refmac(result.structure, cycles=10, auto_accept=False)

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
        print(name, flush=True)
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
        print(json.dumps(result_dict, indent=4), flush=True)
        self.report["jobs"].append({"name": name, **result_dict})
        self.write_report()

    def write_report(self):
        self.seconds["total"] = time.time() - self.start_time
        with open(self.path("modelcraft.json"), "w") as report_file:
            json.dump(self.report, report_file, indent=4)

    def print_refmac_result(self, result: RefmacResult):
        model_stats = ModelStats(result.structure)
        print(f"\nResidues: {model_stats.residues:6d}", flush=True)
        if self.args.mode == "xray":
            print(f"Waters:   {model_stats.waters:6d}", flush=True)
            print(f"R-work:   {result.rwork:6.4f}", flush=True)
            print(f"R-free:   {result.rfree:6.4f}", flush=True)
        if self.args.mode == "em":
            print(f"FSC:      {result.fsc:6.4f}", flush=True)

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
            print("The model cell is incompatible with the data cell", flush=True)
            cell1 = " ".join(f"{x:7.2f}" for x in structure.cell.parameters)
            cell2 = " ".join(f"{x:7.2f}" for x in mtz.cell.parameters)
            print(f"Model: {cell1}  {structure_spacegroup.hm}", flush=True)
            print(f"Data:  {cell2}  {mtz.spacegroup.hm}", flush=True)
            print("Molecular replacement should be used first", flush=True)
            self.terminate("Model cell is incompatible")
        remove_scale(structure=structure)
        if distortion > 0:
            update_cell(structure=structure, new_cell=mtz.cell)
