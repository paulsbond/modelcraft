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
        self.maps = {}
        self.fmean = None
        self.phases = None

    def run(self):
        print(f"# ModelCraft {__version__}", flush=True)
        os.makedirs(self.args.directory, exist_ok=self.args.overwrite_directory)
        self.start_time = time.time()
        self._read_input_maps()
        self._trim_input_maps()
        self._calculate_fmean_and_phases()
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

    def _read_input_maps(self):
        if self.args.half_maps:
            self.maps["half_map1"] = read_map(self.args.half_maps[0])
            self.maps["half_map2"] = read_map(self.args.half_maps[1])
        else:
            self.maps["single_map"] = read_map(self.args.single_map)
        if self.args.build_map:
            self.maps["build_map"] = read_map(self.args.build_map)

    def _trim_input_maps(self):
        if self.args.mask:
            if self.args.mask == "auto":
                map_for_mask = self.maps.get("half_map1", self.maps.get("single_map"))
                mask = EmdaMapMask(map_for_mask).run(self).mask
            else:
                mask = read_map(self.args.mask)
            trimmed = ServalcatTrim(mask, self.maps).run(self)
            self.maps.update(trimmed.maps)

    def _calculate_fmean_and_phases(self):
        if self.args.build_map:
            refmac = RefmacMapToMtz(
                density=self.maps["build_map"],
                resolution=self.args.resolution,
            ).run(self)
            fphi = refmac.fphi
        elif self.args.half_maps:
            nemap = ServalcatNemap(
                halfmap1=self.maps["half_map1"],
                halfmap2=self.maps["half_map2"],
                resolution=self.args.resolution,
            ).run(self)
            fphi = nemap.fphi
        else:
            refmac = RefmacMapToMtz(
                density=self.maps["single_map"],
                resolution=self.args.resolution,
            ).run(self)
            fphi = refmac.fphi
        self.fmean, self.phases = convert_to_fsigf_and_phifom(fphi)

    def buccaneer(self, structure: gemmi.Structure) -> gemmi.Structure:
        result = Buccaneer(
            contents=self.args.contents,
            fsigf=self.fmean,
            phases=self.phases,
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
            fsigf=self.fmean,
            phases=self.phases,
            structure=structure,
        ).run(self)
        return result.structure

    def servalcat_refine(self, structure: gemmi.Structure) -> gemmi.Structure:
        if ModelStats(structure).residues == 0:
            self.terminate(reason="No residues to refine")
        result = ServalcatRefine(
            structure=structure,
            resolution=self.args.resolution,
            halfmap1=self.maps.get("half_map1"),
            halfmap2=self.maps.get("half_map2"),
            density=self.maps.get("single_map"),
            ligand=self.args.restraints,
        ).run(self)
        return result.structure

    def servalcat_fsc(self, structure: gemmi.Structure) -> float:
        result = ServalcatFsc(
            structure=structure,
            resolution=self.args.resolution,
            halfmap1=self.maps.get("half_map1"),
            halfmap2=self.maps.get("half_map2"),
            density=self.maps.get("single_map"),
        ).run(self)
        return result.fsc
