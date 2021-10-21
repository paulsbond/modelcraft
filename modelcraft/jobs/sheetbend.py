import dataclasses
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class SheetbendResult:
    structure: gemmi.Structure
    seconds: float


class Sheetbend(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        fsigf: DataItem,
        freer: DataItem = None,
        regularise: bool = False,
    ):
        super().__init__("csheetbend")
        self.structure = structure
        self.fsigf = fsigf
        self.freer = freer
        self.regularise = regularise

    def _setup(self) -> None:
        write_mmcif(self._path("xyzin.cif"), self.structure)
        write_mtz(self._path("hklin.mtz"), [self.fsigf, self.freer])
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        if self.freer is not None:
            self._args += ["-colin-free", self.freer.label()]
        self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cycles", "12"]
        self._args += ["-resolution-by-cycle", "6,3"]
        if self.regularise:
            self._args += ["-postrefine-u-iso"]
            self._args += ["-pseudo-regularize"]
            self._args += ["-refine-regularize-cycles", "3"]

    def _result(self) -> SheetbendResult:
        self._check_files_exist("xyzout.cif")
        return SheetbendResult(
            structure=read_structure(self._path("xyzout.cif")),
            seconds=self._seconds,
        )
