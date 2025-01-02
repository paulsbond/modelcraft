import dataclasses
import xml.etree.ElementTree as ET
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, remove_non_library_atoms, write_mmcif


@dataclasses.dataclass
class NucleoFindResult:
    predicted_phosphate_map: gemmi.Ccp4Map
    predicted_sugar_map: gemmi.Ccp4Map
    predicted_base_map: gemmi.Ccp4Map
    seconds: float


class NucleoFind(Job):
    def __init__(
        self,
        fphi: DataItem = None
    ):
        super().__init__("nucleofind")
        self.fphi = fphi

    def _setup(self) -> None:
        data_items = [self.fphi]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-i", "hklin.mtz"]
        self._args += ["-o", "nucleofind"]
        self._args += ["-amplitude", self.fphi.label(0)]
        self._args += ["-phase", self.fphi.label(1)]
        self._args += ["-m", "core"]
        self._args += ["-no-symmetry"]
        self._args += ["-gpu"] # If we can know what system they are on, we can use the GPU
        

    def _result(self) -> NucleoFindResult:
        self._check_files_exist("nucleofind/nucleofind-phosphate.map")
        self._check_files_exist("nucleofind/nucleofind-sugar.map")
        self._check_files_exist("nucleofind/nucleofind-base.map")

        predicted_phosphate_map = gemmi.read_ccp4_map(self._path("nucleofind/nucleofind-phosphate.map"))
        predicted_sugar_map = gemmi.read_ccp4_map(self._path("nucleofind/nucleofind-sugar.map"))
        predicted_base_map = gemmi.read_ccp4_map(self._path("nucleofind/nucleofind-base.map"))

        return NucleoFindResult(
            predicted_phosphate_map=predicted_phosphate_map,
            predicted_sugar_map=predicted_sugar_map,
            predicted_base_map=predicted_base_map,
            seconds=self._seconds,
        )
