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
        data_items = [self.fsigf, self.freer, self.phases, self.fphi]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        self._args += ["-colin-free", self.freer.label()]
        if self.phases.types == "AAAA":
            self._args += ["-colin-hl", self.phases.label()]
        else:
            self._args += ["-colin-phifom", self.phases.label()]
        if self.fphi is not None:
            self._args += ["-colin-fc", self.fphi.label()]
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
