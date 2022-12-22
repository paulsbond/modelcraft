import dataclasses
import gemmi
from ..job import Job
from ..maps import write_map
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
        write_map(self._path("map.cpp4"), self.map)
        write_map(self._path("mask.ccp4"), self.mask)
        write_mmcif(self._path("model.cif"), self.structure)
        self._args += ["-m", "servalcat.command_line", "trim"]
        self._args += ["--maps", "map.cpp4"]
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
