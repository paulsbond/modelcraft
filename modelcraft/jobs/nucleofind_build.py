import dataclasses
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, remove_non_library_atoms, write_mmcif


@dataclasses.dataclass
class NucleoFindBuildResult:
    structure: gemmi.Structure
    seconds: float


class NucleoFindBuild(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        freer: DataItem = None,
        structure: gemmi.Structure = None,
        predicted_map: gemmi.Ccp4Map = None,
        cycles: int = 3,
    ):
        super().__init__("nucleofind-build")
        self.contents = contents
        self.fsigf = fsigf
        self.phases = phases
        self.fphi = fphi
        self.freer = freer
        self.structure = structure
        self.cycles = cycles
        self.predicted_map = predicted_map

    def _setup(self) -> None:
        types = [PolymerType.RNA, PolymerType.DNA]
        self.contents.write_sequence_file(self._path("seqin.seq"), types)
        self._args += ["-seqin", "seqin.seq"]
        data_items = [self.fsigf, self.phases, self.fphi, self.freer]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        self.predicted_map.write_ccp4_map(self._path("predin.map"))
        self._args += ["-predin", "predin.map"]
        if self.fphi is not None:
            self._args += ["-colin-fc", self.fphi.label()]
        if self.freer is not None:
            self._args += ["-colin-free", self.freer.label()]
        if self.structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.structure)
            self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-cycles", str(self.cycles)]
        self._args += ["-pdbout", "xyzout.pdb"]

    def _result(self) -> NucleoFindBuildResult:
        self._check_files_exist( "xyzout.pdb")
        structure = read_structure(self._path("xyzout.pdb"))
        return NucleoFindBuildResult(
            structure=structure,
            seconds=self._seconds,
        )
