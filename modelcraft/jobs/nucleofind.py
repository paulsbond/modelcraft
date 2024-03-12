import dataclasses
import xml.etree.ElementTree as ET
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, remove_non_library_atoms, write_mmcif


@dataclasses.dataclass
class NucleoFindResult:
    predicted_map: gemmi.Ccp4Map
    seconds: float


class NucleoFind(Job):
    def __init__(
        self,
        fphi: DataItem = None,
    ):
        super().__init__("nucleofind")
        self.fphi = fphi

    def _setup(self) -> None:
        data_items = [self.fphi]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-i", "hklin.mtz"]
        self._args += ["-o", "pred.map"]
        self._args += ["-intensity", self.fphi.label(0)]
        self._args += ["-phase", self.fphi.label(1)]
        self._args += ["-model_path", "/Users/dialpuri/Development/nucleofind/models/phosphate.onnx"]

    def _result(self) -> NucleoFindResult:
        self._check_files_exist("pred.map")
        predicted_map = gemmi.read_ccp4_map(self._path("pred.map"))
        return NucleoFindResult(
            predicted_map=predicted_map,
            seconds=self._seconds,
        )
