import dataclasses
import shutil
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz


@dataclasses.dataclass
class NucleoFindResult:
    predicted_phosphate_map: gemmi.Ccp4Map
    predicted_sugar_map: gemmi.Ccp4Map
    predicted_base_map: gemmi.Ccp4Map
    seconds: float


class NucleoFind(Job):
    def __init__(
        self,
        fphi: DataItem = None,
        fsigf: DataItem = None,
        phases: DataItem = None,
    ):
        super().__init__("nucleofind")
        self.fphi = fphi
        self.fsigf = fsigf
        self.phases = phases

    @staticmethod
    def is_available():
        if shutil.which("nucleofind") is None:
            raise False
        return True

    def _setup(self) -> None:
        data_items = [self.fphi]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-i", "hklin.mtz"]
        self._args += ["-o", "nucleofind"]
        if self.fphi is not None:
            self._args += ["--amplitude", self.fphi.label(0)]
            self._args += ["--phase", self.fphi.label(1)]
        elif self.fsigf is not None and self.phases is not None:
            self._args += ["--amplitude", self.fsigf.label(0)]
            self._args += ["--phase", self.phases.label(0)]
        else:
            raise RuntimeError("Incorrect columns passed to NucleoFind")
        self._args += ["-m", "core"]
        self._args += ["--gpu"] # If we can know what system they are on, we can use the GPU
        

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
