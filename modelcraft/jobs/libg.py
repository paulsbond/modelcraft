import dataclasses
import gemmi
from ..job import Job
from ..structure import write_mmcif


@dataclasses.dataclass
class LibgResult:
    seconds: float


class Libg(Job):
    def __init__(self, structure: gemmi.Structure):
        super().__init__("libg")
        self.structure = structure

    def _setup(self) -> None:
        write_mmcif(self._path("structure.cif"), self.structure)
        self._args += ["-c", "structure.cif"]
        self._args += ["-o", "libg_restraints.txt"]

    def _result(self) -> LibgResult:
        self._check_files_exist("libg_restraints.txt")
        return LibgResult(
            seconds=self._seconds,
        )
