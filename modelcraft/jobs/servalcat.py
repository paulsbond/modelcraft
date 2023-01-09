import dataclasses
import json
import gemmi
from ..job import Job
from ..maps import write_map
from ..reflections import DataItem
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class ServalcatTrimResult:
    map: gemmi.FloatGrid
    mask: gemmi.FloatGrid
    structure: gemmi.Structure
    seconds: float


class ServalcatTrim(Job):
    def __init__(
        self,
        map_: gemmi.FloatGrid,
        mask: gemmi.FloatGrid,
        structure: gemmi.Structure,
    ):
        super().__init__("ccpem-python")
        self.map = map_
        self.mask = mask
        self.structure = structure

    def _setup(self) -> None:
        write_map(self._path("map.ccp4"), self.map)
        write_map(self._path("mask.ccp4"), self.mask)
        write_mmcif(self._path("model.cif"), self.structure)
        self._args += ["-m", "servalcat.command_line", "trim"]
        self._args += ["--maps", "map.ccp4"]
        self._args += ["--mask", "mask.ccp4"]
        self._args += ["--model", "model.cif"]

    def _result(self) -> ServalcatTrimResult:
        self._check_files_exist(
            "map_trimmed.mrc", "mask_trimmed.mrc", "model_trimmed.cif"
        )
        return ServalcatTrimResult(
            map=gemmi.read_ccp4_map(self._path("map_trimmed.mrc")).grid,
            mask=gemmi.read_ccp4_map(self._path("mask_trimmed.mrc")).grid,
            structure=read_structure(self._path("model_trimmed.cif")),
            seconds=self._seconds,
        )


@dataclasses.dataclass
class ServalcatNemapResult:
    fphi: DataItem
    seconds: float


class ServalcatNemap(Job):
    def __init__(
        self,
        halfmap1: gemmi.Ccp4Map,
        halfmap2: gemmi.Ccp4Map,
        mask: gemmi.FloatGrid,
        resolution: float,
    ):
        super().__init__("ccpem-python")
        self.halfmap1 = halfmap1
        self.halfmap2 = halfmap2
        self.mask = mask
        self.resolution = resolution

    def _setup(self) -> None:
        self.halfmap1.write_ccp4_map(self._path("halfmap1.ccp4"))
        self.halfmap2.write_ccp4_map(self._path("halfmap2.ccp4"))
        write_map(self._path("mask.ccp4"), self.mask)
        self._args += ["-m", "servalcat.command_line", "util", "nemap"]
        self._args += ["--halfmaps", "halfmap1.ccp4", "halfmap2.ccp4"]
        self._args += ["--mask", "mask.ccp4"]
        self._args += ["--resolution", str(self.resolution)]

    def _result(self) -> ServalcatNemapResult:
        self._check_files_exist("nemap.mtz")
        mtz = gemmi.read_mtz_file(self._path("nemap.mtz"))
        return ServalcatNemapResult(
            fphi=DataItem(mtz, "FWT,PHWT"),
            seconds=self._seconds,
        )


@dataclasses.dataclass
class ServalcatRefineResult:
    structure: gemmi.Structure
    fphi_best: DataItem
    fphi_diff: DataItem
    fphi_calc: DataItem
    fsc: float
    seconds: float


class ServalcatRefine(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        density: gemmi.FloatGrid,
        resolution: float,
        blur: float = 0.0,
        cycles: int = 20,
        bfactor: float = 40.0,
    ):
        super().__init__("ccpem-python")
        self.structure = structure
        self.density = density
        self.resolution = resolution
        self.blur = blur
        self.cycles = cycles
        self.bfactor = bfactor

    def _setup(self) -> None:
        write_mmcif(self._path("model.cif"), self.structure)
        write_map(self._path("map.ccp4"), self.density)
        self._args += ["-m", "servalcat.command_line", "refine_spa"]
        self._args += ["--model", "model.cif"]
        self._args += ["--map", "map.ccp4"]
        self._args += ["--no_mask"]
        self._args += ["--blur", str(self.blur)]
        self._args += ["--ncycle", str(self.cycles)]
        self._args += ["--bfactor", str(self.bfactor)]
        self._args += ["--resolution", str(self.resolution)]

    def _result(self) -> ServalcatRefineResult:
        self._check_files_exist("refined.mmcif", "refined.mtz", "refined_summary.json")
        mtz = gemmi.read_mtz_file(self._path("refined.mtz"))
        with open(self._path("refined_summary.json")) as json_file:
            summary = json.load(json_file)
        return ServalcatRefineResult(
            structure=read_structure(self._path("refined.mmcif")),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fphi_calc=DataItem(mtz, "FC_ALL,PHIC_ALL"),
            fsc=float(summary["cycles"][-1]["fsc_average"]),
            seconds=self._seconds,
        )
