import dataclasses
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class NautilusResult:
    structure: gemmi.Structure


class Nautilus(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        structure: gemmi.Structure = None,
        cycles: int = 3,
    ):
        super().__init__("cnautilus")
        self.contents = contents
        self.fsigf = fsigf
        self.freer = freer
        self.phases = phases
        self.fphi = fphi
        self.structure = structure
        self.cycles = cycles

    def _setup(self) -> None:
        types = [PolymerType.RNA, PolymerType.DNA]
        self.contents.write_sequence_file(self._path("seqin.seq"), types)
        self._args += ["-seqin", "seqin.seq"]
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
            self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-cycles", str(self.cycles)]
        self._args += ["-anisotropy-correction"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cif"]

    def _result(self) -> NautilusResult:
        return NautilusResult(structure=read_structure(self._path("xyzout.cif")))
