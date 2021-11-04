import dataclasses
import gemmi
from ..job import Job
from ..reflections import DataItem


@dataclasses.dataclass
class FreeRFlagResult:
    freer: DataItem
    seconds: float


class FreeRFlag(Job):
    def __init__(self, mtz: gemmi.Mtz):
        super().__init__("freerflag")
        self.mtz = mtz

    def _setup(self) -> None:
        self.mtz.write_to_file(self._path("hklin.mtz"))
        self._args += ["HKLIN", "hklin.mtz"]
        self._args += ["HKLOUT", "hklout.mtz"]
        self._stdin.append("UNIQUE")
        self._stdin.append("END")

    def _result(self) -> FreeRFlagResult:
        self._check_files_exist("hklout.mtz")
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        freer = DataItem(mtz, "FreeR_flag")
        return FreeRFlagResult(
            freer=freer,
            seconds=self._seconds,
        )
