import dataclasses
import xml.etree.ElementTree as ET
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, remove_non_library_atoms, write_mmcif


@dataclasses.dataclass
class NautilusResult:
    structure: gemmi.Structure
    fragments_built: int
    residues_built: int
    residues_sequenced: int
    longest_fragment: int
    seconds: float


class Nautilus(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        freer: DataItem = None,
        structure: gemmi.Structure = None,
        cycles: int = 3,
    ):
        super().__init__("cnautilus")
        self.contents = contents
        self.fsigf = fsigf
        self.phases = phases
        self.fphi = fphi
        self.freer = freer
        self.structure = structure
        self.cycles = cycles

    def _setup(self) -> None:
        types = [PolymerType.RNA, PolymerType.DNA]
        self.contents.write_sequence_file(self._path("seqin.seq"), types)
        self._args += ["-seqin", "seqin.seq"]
        data_items = [self.fsigf, self.phases, self.fphi, self.freer]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        if self.phases.types == "AAAA":
            self._args += ["-colin-hl", self.phases.label()]
        else:
            self._args += ["-colin-phifom", self.phases.label()]
        if self.fphi is not None:
            self._args += ["-colin-fc", self.fphi.label()]
        if self.freer is not None:
            self._args += ["-colin-free", self.freer.label()]
        if self.structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.structure)
            self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-cycles", str(self.cycles)]
        self._args += ["-anisotropy-correction"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cif"]
        self._args += ["-xmlout", "xmlout.xml"]

    def _result(self) -> NautilusResult:
        self._check_files_exist("xmlout.xml", "xyzout.cif")
        xml = ET.parse(self._path("xmlout.xml")).getroot()
        structure = read_structure(self._path("xyzout.cif"))
        remove_non_library_atoms(structure)
        return NautilusResult(
            structure=structure,
            fragments_built=int(xml.find("Final/FragmentsBuilt").text),
            residues_built=int(xml.find("Final/ResiduesBuilt").text),
            residues_sequenced=int(xml.find("Final/ResiduesSequenced").text),
            longest_fragment=int(xml.find("Final/ResiduesLongestFragment").text),
            seconds=self._seconds,
        )
