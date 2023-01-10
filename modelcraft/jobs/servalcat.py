import dataclasses
import json
import gemmi
from ..job import Job
from ..reflections import DataItem
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class ServalcatNemapResult:
    fphi: DataItem
    seconds: float


class ServalcatNemap(Job):
    def __init__(
        self,
        halfmap1: gemmi.Ccp4Map,
        halfmap2: gemmi.Ccp4Map,
        mask: gemmi.Ccp4Map,
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
        self.mask.write_ccp4_map(self._path("mask.ccp4"))
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
class ServalcatTrimResult:
    mask: gemmi.Ccp4Map
    density: gemmi.Ccp4Map
    halfmap1: gemmi.Ccp4Map
    halfmap2: gemmi.Ccp4Map
    structure: gemmi.Structure
    seconds: float


class ServalcatTrim(Job):
    def __init__(
        self,
        mask: gemmi.Ccp4Map,
        density: gemmi.Ccp4Map = None,
        halfmap1: gemmi.Ccp4Map = None,
        halfmap2: gemmi.Ccp4Map = None,
        structure: gemmi.Structure = None,
    ):
        super().__init__("ccpem-python")
        self.mask = mask
        self.density = density
        self.halfmap1 = halfmap1
        self.halfmap2 = halfmap2
        self.structure = structure

    def _setup(self) -> None:
        self._args += ["-m", "servalcat.command_line", "trim"]
        self.mask.write_ccp4_map(self._path("mask.ccp4"))
        self._args += ["--mask", "mask.ccp4"]
        if any(x is not None for x in (self.density, self.halfmap1, self.halfmap2)):
            self._args.append("--maps")
        if self.density is not None:
            self.density.write_ccp4_map(self._path("density.ccp4"))
            self._args.append("density.ccp4")
        if self.halfmap1 is not None:
            self.halfmap1.write_ccp4_map(self._path("halfmap1.ccp4"))
            self._args.append("halfmap1.ccp4")
        if self.halfmap2 is not None:
            self.halfmap2.write_ccp4_map(self._path("halfmap2.ccp4"))
            self._args.append("halfmap2.ccp4")
        if self.structure is not None:
            write_mmcif(self._path("model.cif"), self.structure)
            self._args += ["--model", "model.cif"]

    def _result(self) -> ServalcatTrimResult:
        return ServalcatTrimResult(
            mask=gemmi.read_ccp4_map(self._path("mask_trimmed.mrc")),
            density=None
            if self.density is None
            else gemmi.read_ccp4_map(self._path("density_trimmed.mrc")),
            halfmap1=None
            if self.halfmap1 is None
            else gemmi.read_ccp4_map(self._path("halfmap1_trimmed.mrc")),
            halfmap2=None
            if self.halfmap2 is None
            else gemmi.read_ccp4_map(self._path("halfmap2_trimmed.mrc")),
            structure=None
            if self.structure is None
            else read_structure(self._path("model_trimmed.cif")),
            seconds=self._seconds,
        )


@dataclasses.dataclass
class ServalcatRefineResult:
    structure: gemmi.Structure
    fphi_best: DataItem
    fphi_diff: DataItem
    fsc: float
    seconds: float


class ServalcatRefine(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        resolution: float,
        halfmap1: gemmi.Ccp4Map = None,
        halfmap2: gemmi.Ccp4Map = None,
        density: gemmi.Ccp4Map = None,
        blur: float = 0.0,
        cycles: int = 20,
        bfactor: float = 40.0,
        jellybody: bool = True,
        jellybody_sigma: float = 0.01,
        jellybody_dmax: float = 4.2,
    ):
        super().__init__("ccpem-python")
        self.structure = structure
        self.resolution = resolution
        self.halfmap1 = halfmap1
        self.halfmap2 = halfmap2
        self.density = density
        self.blur = blur
        self.cycles = cycles
        self.bfactor = bfactor
        self.jellybody = jellybody
        self.jellybody_sigma = jellybody_sigma
        self.jellybody_dmax = jellybody_dmax

    def _setup(self) -> None:
        self._args += ["-m", "servalcat.command_line", "refine_spa"]
        write_mmcif(self._path("structure.cif"), self.structure)
        self._args += ["--model", "structure.cif"]
        if self.halfmap1 is not None and self.halfmap2 is not None:
            self.halfmap1.write_ccp4_map(self._path("halfmap1.ccp4"))
            self.halfmap2.write_ccp4_map(self._path("halfmap2.ccp4"))
            self._args += ["--halfmaps", "halfmap1.ccp4", "halfmap2.ccp4"]
        else:
            self.density.write_ccp4_map(self._path("density.ccp4"))
            self._args += ["--map", "density.ccp4"]
        self._args += ["--resolution", str(self.resolution)]
        self._args += ["--blur", str(self.blur)]
        self._args += ["--ncycle", str(self.cycles)]
        self._args += ["--bfactor", str(self.bfactor)]
        self._args += ["--no_mask"]
        if self.jellybody:
            self._args += [
                "--jellybody",
                "--jellybody_params",
                str(self.jellybody_sigma),
                str(self.jellybody_dmax),
            ]

    def _result(self) -> ServalcatRefineResult:
        self._check_files_exist("refined.mmcif", "refined.mtz", "refined_summary.json")
        mtz = gemmi.read_mtz_file(self._path("refined.mtz"))
        with open(self._path("refined_summary.json")) as json_file:
            summary = json.load(json_file)
        return ServalcatRefineResult(
            structure=read_structure(self._path("refined.mmcif")),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fsc=float(summary["cycles"][-1]["fsc_average"]),
            seconds=self._seconds,
        )
