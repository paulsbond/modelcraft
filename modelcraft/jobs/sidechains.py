import dataclasses
import os

import gemmi

from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class SidechainsResult:
    structure: gemmi.Structure
    seconds: float


class Sidechains(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem):
        super().__init__("modelcraft-sidechains")
        self.structure = structure
        self.fphi = fphi

    def _setup(self) -> None:
        write_mmcif(self._path("xyzin.cif"), self.structure)
        write_mtz(self._path("hklin.mtz"), [self.fphi], ["FWT,PHWT"])
        self._args += ["xyzin.cif", "hklin.mtz", "xyzout.cif"]

    def _result(self) -> SidechainsResult:
        path = self._path("xyzout.cif")
        structure = read_structure(path) if os.path.exists(path) else None
        return SidechainsResult(
            structure=structure,
            seconds=self._seconds,
        )
