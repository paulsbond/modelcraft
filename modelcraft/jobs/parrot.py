import dataclasses
import gemmi
from ..contents import AsuContents
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..solvent import solvent_fraction
from ..structure import write_mmcif


@dataclasses.dataclass
class ParrotResult:
    abcd: DataItem
    fphi: DataItem
    seconds: float


class Parrot(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        structure: gemmi.Structure = None,
    ):
        super().__init__("cparrot")
        self.contents = contents
        self.fsigf = fsigf
        self.freer = freer
        self.phases = phases
        self.fphi = fphi
        self.structure = structure

    def _setup(self) -> None:
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", "F,SIGF"]
        self._args += ["-colin-free", "FREE"]
        phases_arg = "-colin-hl" if self.phases.types == "AAAA" else "-colin-phifom"
        phases_label = "HLA,HLB,HLC,HLD" if self.phases.types == "AAAA" else "PHI,FOM"
        self._args += [phases_arg, phases_label]
        if self.fphi is not None:
            self._args += ["-colin-fc", "FC,PHIC"]
        write_mtz(
            path=self._path("hklin.mtz"),
            items=[self.fsigf, self.freer, self.phases, self.fphi],
            labels=["F,SIGF", "FREE", phases_label, "FC,PHIC"],
        )
        if self.structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.structure)
            self._args += ["-pdbin-mr", "xyzin.cif"]
        solvent_content = solvent_fraction(self.contents, self.fsigf)
        self._args += ["-solvent-content", "%.3f" % solvent_content]
        self._args += ["-cycles", "5"]
        self._args += ["-anisotropy-correction"]
        self._args += ["-mtzout", "hklout.mtz"]

    def _result(self) -> ParrotResult:
        self._check_files_exist("hklout.mtz")
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        return ParrotResult(
            abcd=DataItem(mtz, "parrot.ABCD"),
            fphi=DataItem(mtz, "parrot.F_phi"),
            seconds=self._seconds,
        )
