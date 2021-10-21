import dataclasses
import xml.etree.ElementTree as ET
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure


@dataclasses.dataclass
class MolrepResult:
    structure: gemmi.Structure
    err_level: int
    err_message: str
    n_solution: int
    mr_score: float
    mr_zscore: float
    seconds: float


class Molrep(Job):
    def __init__(
        self,
        observations: DataItem,
        structure: gemmi.Structure,
        number_of_monomers: int = 1,
    ):
        super().__init__("molrep")
        self.observations = observations
        self.structure = structure
        self.number_of_monomers = number_of_monomers

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.observations])
        self.structure.write_minimal_pdb(self._path("xyzin.pdb"))
        self._args += ["-f", "hklin.mtz"]
        self._args += ["-m", "xyzin.pdb"]
        self._args += ["-i"]
        labin = "LABIN"
        if self.observations.types == "FQ":
            labin += " F=" + self.observations.label(0)
            labin += " SIGF=" + self.observations.label(1)
        elif self.observations.types == "GLGL":
            labin += " F=" + self.observations.label(0)
            labin += " SIGF=" + self.observations.label(1)
            labin += " F(-)=" + self.observations.label(2)
            labin += " SIGF(-)=" + self.observations.label(3)
        elif self.observations.types == "JQ":
            labin += " I=" + self.observations.label(0)
            labin += " SIGI=" + self.observations.label(1)
        elif self.observations.types == "KMKM":
            labin += " I=" + self.observations.label(0)
            labin += " SIGI=" + self.observations.label(1)
            labin += " I(-)=" + self.observations.label(2)
            labin += " SIGI(-)=" + self.observations.label(3)
        else:
            message = "Unknown MTZ column types for observations: "
            message += self.observations.types
            raise RuntimeError(message)
        self._stdin += [labin]
        self._stdin += [f"NMON {self.number_of_monomers}"]

    def _result(self) -> MolrepResult:
        self._check_files_exist("molrep.pdb", "molrep.xml")
        xml = ET.parse(self._path("molrep.xml")).getroot()
        return MolrepResult(
            structure=read_structure(self._path("molrep.pdb")),
            err_level=int(xml.find("err_level").text),
            err_message=xml.find("err_message").text.strip(),
            n_solution=int(xml.find("n_solution").text),
            mr_score=float(xml.find("mr_score").text),
            mr_zscore=float(xml.find("mr_zscore").text),
            seconds=self._seconds,
        )
