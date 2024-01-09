import os
import time
import gemmi
from . import __version__
from .jobs.buccaneer import Buccaneer
from .jobs.emda import EmdaMapMask
from .jobs.nautilus import Nautilus
from .jobs.refmac import RefmacMapToMtz
from .jobs.servalcat import ServalcatFsc, ServalcatNemap, ServalcatRefine, ServalcatTrim
from .maps import read_map
from .pipeline import Pipeline
from .reflections import convert_to_fsigf_and_phifom
from .structure import ModelStats, write_mmcif


class ModelCraftEm(Pipeline):
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

    def run(self):
        print(f"# ModelCraft {__version__}", flush=True)
        os.makedirs(self.args.directory, exist_ok=self.args.overwrite_directory)
        self.start_time = time.time()
        self._process_input_maps()
        structure = self.args.model
        best_fsc = None
        cycles_without_improvement = 0
        for cycle in range(1, self.args.cycles + 1):
            print(f"\n## Cycle {cycle}\n", flush=True)
            if self.args.contents.proteins:
                structure = self.buccaneer(structure)
                structure = self.servalcat_refine(structure)
            if self.args.contents.rnas or self.args.contents.dnas:
                structure = self.nautilus(structure)
                structure = self.servalcat_refine(structure)
            model_stats = ModelStats(structure)
            fsc = self.servalcat_fsc(structure)
            stats = {"cycle": cycle, "residues": model_stats.residues, "fsc": fsc}
            self.report["cycles"].append(stats)
            if best_fsc is None or fsc > best_fsc:
                best_fsc = fsc
                cycles_without_improvement = 0
                write_mmcif(self.path("modelcraft.cif"), structure)
                self.report["final"] = stats
            else:
                cycles_without_improvement += 1
            self.write_report()
            if cycles_without_improvement == self.args.auto_stop_cycles > 0:
                break
        self.terminate("Normal")

    def _process_input_maps(self):
        maps = [read_map(path) for path in self.args.map]
        if self.args.mask is not None:
            mask = read_map(self.args.mask)
        else:
            mask = EmdaMapMask(maps[0]).run(self).mask
        trimmed = ServalcatTrim(mask, maps).run(self)
        self.args.map = trimmed.maps
        if len(maps) == 2:
            nemap = ServalcatNemap(
                halfmap1=trimmed.maps[0],
                halfmap2=trimmed.maps[1],
                resolution=self.args.resolution,
                mask=trimmed.mask,
            ).run(self)
            self.args.fphi = nemap.fphi
        else:
            refmac = RefmacMapToMtz(
                density=trimmed.maps[0],
                resolution=self.args.resolution,
                blur=self.args.blur,
            ).run(self)
            self.args.fphi = refmac.fphi
        self.args.fmean, self.args.phases = convert_to_fsigf_and_phifom(self.args.fphi)

    def buccaneer(self, structure: gemmi.Structure) -> gemmi.Structure:
        result = Buccaneer(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.args.phases,
            input_structure=structure,
            mr_structure=self.args.model,
            cycles=5,
            threads=self.args.threads,
            em_mode=True,
        ).run(self)
        return result.structure

    def nautilus(self, structure: gemmi.Structure) -> gemmi.Structure:
        result = Nautilus(
            contents=self.args.contents,
            fsigf=self.args.fmean,
            phases=self.args.phases,
            structure=structure,
        ).run(self)
        return result.structure

    def servalcat_refine(self, structure: gemmi.Structure) -> gemmi.Structure:
        if ModelStats(structure).residues == 0:
            self.terminate(reason="No residues to refine")
        result = ServalcatRefine(
            structure=structure,
            resolution=self.args.resolution,
            halfmap1=self.args.map[0] if len(self.args.map) == 2 else None,
            halfmap2=self.args.map[1] if len(self.args.map) == 2 else None,
            density=self.args.map[0] if len(self.args.map) == 1 else None,
            ligand=self.args.restraints,
        ).run(self)
        return result.structure

    def servalcat_fsc(self, structure: gemmi.Structure) -> float:
        result = ServalcatFsc(
            structure=structure,
            resolution=self.args.resolution,
            halfmap1=self.args.map[0] if len(self.args.map) == 2 else None,
            halfmap2=self.args.map[1] if len(self.args.map) == 2 else None,
            density=self.args.map[0] if len(self.args.map) == 1 else None,
        ).run(self)
        return result.fsc
