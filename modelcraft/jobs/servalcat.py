import dataclasses
import shutil
import gemmi
from ..job import Job
from ..maps import read_map
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
        resolution: float,
        mask: gemmi.Ccp4Map = None,
    ):
        super().__init__("servalcat")
        self.halfmap1 = halfmap1
        self.halfmap2 = halfmap2
        self.mask = mask
        self.resolution = resolution

    def _setup(self) -> None:
        self.halfmap1.write_ccp4_map(self._path("halfmap1.ccp4"))
        self.halfmap2.write_ccp4_map(self._path("halfmap2.ccp4"))
        self._args += ["util", "nemap"]
        self._args += ["--halfmaps", "halfmap1.ccp4", "halfmap2.ccp4"]
        self._args += ["--resolution", str(self.resolution)]
        if self.mask is not None:
            self.mask.write_ccp4_map(self._path("mask.ccp4"))
            self._args += ["--mask", "mask.ccp4"]

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
    maps: list
    seconds: float


class ServalcatTrim(Job):
    def __init__(self, mask: gemmi.Ccp4Map, maps: list):
        super().__init__("servalcat")
        self.mask = mask
        self.maps = maps

    def _setup(self) -> None:
        self._args += ["trim"]
        self.mask.write_ccp4_map(self._path("mask.ccp4"))
        self._args += ["--mask", "mask.ccp4"]
        self._args.append("--maps")
        for i, map_ in enumerate(self.maps):
            map_.write_ccp4_map(self._path(f"map{i}.ccp4"))
            self._args.append(f"map{i}.ccp4")
        self._args.append("--noncubic")
        self._args.append("--noncentered")
        self._args.append("--no_shift")

    def _result(self) -> ServalcatTrimResult:
        self._check_files_exist("mask_trimmed.mrc", "map0_trimmed.mrc")
        return ServalcatTrimResult(
            mask=read_map(self._path("mask_trimmed.mrc")),
            maps=[
                read_map(self._path(f"map{i}_trimmed.mrc"))
                for i in range(len(self.maps))
            ],
            seconds=self._seconds,
        )


@dataclasses.dataclass
class ServalcatRefineResult:
    structure: gemmi.Structure
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
        ligand: str = None,
    ):
        super().__init__("servalcat")
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
        self.ligand = ligand

    def _setup(self) -> None:
        self._args += ["refine_spa"]
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
        if self.jellybody:
            self._args += [
                "--jellybody",
                "--jellybody_params",
                str(self.jellybody_sigma),
                str(self.jellybody_dmax),
            ]
        if self.ligand:
            shutil.copy(self.ligand, self._path("ligand.cif"))
            self._args += ["--ligand", "ligand.cif"]

    def _result(self) -> ServalcatRefineResult:
        self._check_files_exist("refined.mmcif")
        return ServalcatRefineResult(
            structure=read_structure(self._path("refined.mmcif")),
            seconds=self._seconds,
        )


@dataclasses.dataclass
class ServalcatFscResult:
    fsc: float
    seconds: float


class ServalcatFsc(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        resolution: float,
        halfmap1: gemmi.Ccp4Map = None,
        halfmap2: gemmi.Ccp4Map = None,
        density: gemmi.Ccp4Map = None,
    ):
        super().__init__("servalcat")
        self.structure = structure
        self.resolution = resolution
        self.halfmap1 = halfmap1
        self.halfmap2 = halfmap2
        self.density = density

    def _setup(self) -> None:
        self._args += ["fsc"]
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

    def _result(self) -> ServalcatFscResult:
        self._check_files_exist("fsc.dat")
        fsc = None
        with open(self._path("fsc.dat")) as text:
            for line in text:
                if line.startswith("# FSCaverage of fsc_FC_full ="):
                    fsc = float(line.strip().split()[-1])
        return ServalcatFscResult(
            fsc=fsc,
            seconds=self._seconds,
        )
