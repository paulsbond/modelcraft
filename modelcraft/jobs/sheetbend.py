import dataclasses
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class SheetbendResult:
    structure: gemmi.Structure


class Sheetbend(Job):
    def __init__(
        self,
        fsigf: DataItem,
        freer: DataItem,
        structure: gemmi.Structure,
        executable: str = None,
    ):
        super().__init__(executable or "csheetbend")
        self.fsigf = fsigf
        self.freer = freer
        self.structure = structure

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fsigf, self.freer])
        write_mmcif(self._path("xyzin.cif"), self.structure)
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        self._args += ["-colin-free", self.freer.label()]
        self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cycles", "12"]
        self._args += ["-resolution-by-cycle", "6,3"]

    def _result(self) -> SheetbendResult:
        return SheetbendResult(structure=read_structure(self._path("xyzout.cif")))
