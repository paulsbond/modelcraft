import dataclasses
import gemmi
from ..job import Job


@dataclasses.dataclass
class AcedrgResult:
    chemcomp: gemmi.ChemComp
    seconds: float


class Acedrg(Job):
    def __init__(self, code: str, cif: gemmi.cif.Document = None, smiles: str = None):
        super().__init__("acedrg")
        self.code = code
        self.cif = cif
        self.smiles = smiles

    def _setup(self) -> None:
        self._args += ["-r", self.code.upper()]
        self._args += ["-o", "output"]
        if self.cif is not None:
            self.cif.write_file(self._path("input.cif"))
            self._args += ["-c", "input.cif"]
        if self.smiles is not None:
            self._args += ["-i", self.smiles]

    def _result(self) -> AcedrgResult:
        self._check_files_exist("output.cif")
        block = gemmi.cif.read(self._path("output.cif"))[-1]
        return AcedrgResult(
            chemcomp=gemmi.make_chemcomp_from_block(block),
            seconds=self._seconds,
        )
