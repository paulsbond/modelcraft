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
        self._args += ["-o", "pred"]
        self._args += ["-intensity", self.fphi.label(0)]
        self._args += ["-phase", self.fphi.label(1)]
        self._args += ["-m", "all"]
        # self._args += ["-gpu"]

    def _result(self) -> NucleoFindResult:
        self._check_files_exist("pred_phosphate.map", "pred_sugar.map", "pred_base.map")

        predicted_phosphate_map = gemmi.read_ccp4_map(self._path("pred_phosphate.map"))
        predicted_sugar_map = gemmi.read_ccp4_map(self._path("pred_sugar.map"))
        predicted_base_map = gemmi.read_ccp4_map(self._path("pred_sugar.map"))

        return NucleoFindResult(
            predicted_phosphate_map=predicted_phosphate_map,
            predicted_sugar_map=predicted_sugar_map,
            predicted_base_map=predicted_base_map,
            seconds=self._seconds,
        )
