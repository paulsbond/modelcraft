import os
import time
import gemmi
from . import __version__
from .jobs.buccaneer import Buccaneer
from .jobs.coot import FixSideChains, Prune
from .jobs.ctruncate import CTruncate
from .jobs.findwaters import FindWaters
from .jobs.nautilus import Nautilus
from .jobs.parrot import Parrot
from .jobs.refmac import Refmac
from .jobs.sheetbend import Sheetbend
from .cell import max_distortion, remove_scale, update_cell
from .pipeline import Pipeline
from .reflections import DataItem, write_mtz
from .structure import ModelStats, remove_residues, write_mmcif


class ModelCraftXray(Pipeline):
    def __init__(self, parsed_args, raw_args):
        self.args = parsed_args
        super().__init__(
            directory=self.args.directory,
            keep_jobs=self.args.keep_files,
            keep_logs=self.args.keep_logs,
            json_name="modelcraft.json",
        )
        self.report["version"] = __version__
        self.report["args"] = raw_args
        self.report["cycles"] = []
        self.cycle = 0
        self.current_structure: gemmi.Structure = self.args.model
        self.current_phases: DataItem = self.args.phases
        self.current_fphi_best: DataItem = None
        self.current_fphi_diff: DataItem = None
        self.current_fphi_calc: DataItem = None
        self.last_refmac = None
        self.output_refmac = None
        self.cycles_without_improvement = 0

    @property
    def resolution(self):
        return self.args.fmean.resolution_high()

    def run(self):
        print(f"# ModelCraft {__version__}", flush=True)
        os.makedirs(self.args.directory, exist_ok=self.args.overwrite_directory)
        self.start_time = time.time()
        if self.args.fmean is None:
            self._convert_observations()
        if self.args.model is not None:
            self._refine_input_model()
        for self.cycle in range(1, self.args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle, flush=True)
            self.run_cycle()
            self.process_cycle_output(self.last_refmac)
            if self.cycles_without_improvement == self.args.auto_stop_cycles > 0:
                break
        if (
            not self.args.basic
            and self.output_refmac.rwork < 0.3
            and self.resolution < 2.5
        ):
            print("\n## Finalisations\n", flush=True)
            self.cycle += 1
            self.update_current_from_refmac_result(self.output_refmac)
            self.fixsidechains()
            self.process_cycle_output(self.last_refmac)
        print("\n## Best Model:", flush=True)
        _print_refmac_result(self.output_refmac)
        self._remove_current_files()
        self.terminate(reason="Normal")

    def _convert_observations(self):
        print("\n## Converting input observations to mean amplitudes\n", flush=True)
        observations = self.args.ianom or self.args.imean or self.args.fanom
        ctruncate = CTruncate(observations=observations).run(self)
        self.args.fmean = ctruncate.fmean
        if self.args.fanom is None and ctruncate.fanom is not None:
            self.args.fanom = ctruncate.fanom
        if self.args.imean is None and ctruncate.imean is not None:
            self.args.imean = ctruncate.imean

    def _refine_input_model(self):
        print("\n## Refining Input Model\n", flush=True)
        self.update_model_cell()
        write_mmcif(self.path("current.cif"), self.current_structure)
        if self.args.disable_sheetbend:
            self.refmac(self.current_structure, cycles=10, auto_accept=True)
        else:
            sheetbend = Sheetbend(
                fsigf=self.args.fmean,
                freer=self.args.freer,
                structure=self.current_structure,
            ).run(self)
            result1 = self.run_refmac(sheetbend.structure, cycles=10)
            result2 = self.run_refmac(self.current_structure, cycles=10)
            chosen = result1 if result1.rfree < result2.rfree else result2
            self.update_current_from_refmac_result(chosen)
        self.args.model = self.current_structure
        if self.args.phases is not None:
            self.current_phases = self.args.phases
        _print_refmac_result(self.last_refmac)

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
                if self.cycle > 1 or self.args.phases is None:
                    self.findwaters(dummy=True)
                remove_residues(structure=self.current_structure, names={"HOH", "DUM"})
            self.buccaneer()
            self.prune(chains_only=True)
            self.nautilus()
            self.findwaters()

    def buccaneer(self):
        if not self.args.contents.proteins:
            return
        result = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            freer=self.args.freer,
            input_structure=self.current_structure,
            mr_structure=self.args.model,
            use_mr=True,
            filter_mr=True,
            seed_mr=True,
            cycles=3 if self.cycle == 1 else 2,
            threads=self.args.threads,
        ).run(self)
        write_mmcif(self.path("current.cif"), result.structure)
        self.refmac(result.structure, cycles=10, auto_accept=True)

    def nautilus(self):
        if not (self.args.contents.rnas or self.args.contents.dnas):
            return
        result = Nautilus(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            freer=self.args.freer,
            structure=self.current_structure,
        ).run(self)
        write_mmcif(self.path("current.cif"), result.structure)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def refmac(self, structure: gemmi.Structure, cycles: int, auto_accept: bool):
        result = self.run_refmac(structure, cycles)
        if auto_accept or result.rfree < self.last_refmac.rfree:
            if not auto_accept:
                print("(accepted)", flush=True)
            self.update_current_from_refmac_result(result)
        else:
            print("(rejected)", flush=True)
            write_mmcif(self.path("current.cif"), self.current_structure)

    def run_refmac(self, structure: gemmi.Structure, cycles: int):
        if ModelStats(structure).residues == 0:
            self.terminate(reason="No residues to refine")
        use_phases = self.args.unbiased and (
            self.output_refmac is None or self.output_refmac.rwork > 0.35
        )
        return Refmac(
            structure=structure,
            fsigf=self.args.fmean,
            freer=self.args.freer,
            cycles=cycles,
            phases=self.args.phases if use_phases else None,
            twinned=self.args.twinned,
            libin=self.args.restraints,
        ).run(self)

    def update_current_from_refmac_result(self, result):
        self.current_structure = result.structure
        self.current_phases = getattr(result, "abcd", None)
        self.current_fphi_best = result.fphi_best
        self.current_fphi_diff = result.fphi_diff
        self.current_fphi_calc = result.fphi_calc
        self.last_refmac = result
        write_mmcif(self.path("current.cif"), result.structure)
        write_mtz(self.path("current.mtz"), [self.current_fphi_best], ["F,PHI"])

    def parrot(self):
        if self.args.disable_parrot:
            return
        result = Parrot(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            freer=self.args.freer,
            phases=self.current_phases,
            fphi=self.current_fphi_best,
            structure=self.current_structure,
        ).run(self)
        self.current_phases = result.abcd
        self.current_fphi_best = result.fphi
        write_mtz(self.path("current.mtz"), [self.current_fphi_best], ["F,PHI"])

    def prune(self, chains_only=False):
        if self.args.disable_pruning or not self.args.contents.proteins:
            return
        result = Prune(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
            chains_only=chains_only,
        ).run(self)
        write_mmcif(self.path("current.cif"), result.structure)
        self.refmac(result.structure, cycles=5, auto_accept=True)

    def fixsidechains(self):
        if self.args.disable_side_chain_fixing or not self.args.contents.proteins:
            return
        result = FixSideChains(
            structure=self.current_structure,
            fphi_best=self.current_fphi_best,
            fphi_diff=self.current_fphi_diff,
        ).run(self)
        write_mmcif(self.path("current.cif"), result.structure)
        self.refmac(result.structure, cycles=5, auto_accept=False)

    def findwaters(self, dummy=False):
        if dummy and self.args.disable_dummy_atoms:
            return
        if not dummy and self.args.disable_waters:
            return
        result = FindWaters(
            structure=self.current_structure,
            fphi=self.current_fphi_best,
            dummy=dummy,
        ).run(self)
        write_mmcif(self.path("current.cif"), result.structure)
        self.refmac(result.structure, cycles=10, auto_accept=False)

    def process_cycle_output(self, result):
        _print_refmac_result(result)
        model_stats = ModelStats(result.structure)
        stats = {"cycle": self.cycle, "residues": model_stats.residues}
        stats["waters"] = model_stats.waters
        stats["r_work"] = result.rwork
        stats["r_free"] = result.rfree
        self.report["cycles"].append(stats)
        if self.output_refmac is None or result.rwork < self.output_refmac.rwork:
            self.cycles_without_improvement = 0
            self.output_refmac = result
            write_mmcif(self.path("modelcraft.cif"), result.structure)
            result.mtz.write_to_file(self.path("modelcraft.mtz"))
            self.report["final"] = stats
        else:
            self.cycles_without_improvement += 1
        self.write_report()

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

    def _remove_current_files(self):
        for filename in ("current.cif", "current.mtz"):
            try:
                os.remove(self.path(filename))
            except FileNotFoundError:
                pass


def _print_refmac_result(result):
    model_stats = ModelStats(result.structure)
    print(f"\nResidues: {model_stats.residues:6d}", flush=True)
    print(f"Waters:   {model_stats.waters:6d}", flush=True)
    print(f"R-work:   {result.rwork:6.4f}", flush=True)
    print(f"R-free:   {result.rfree:6.4f}", flush=True)
