import dataclasses
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz


@dataclasses.dataclass
class ComitResult:
    abcd: DataItem
    fphi: DataItem
    seconds: float


class Comit(Job):
    def __init__(self, fsigf: DataItem, fphi: DataItem):
        super().__init__("comit")
        self.fsigf = fsigf
        self.fphi = fphi

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fsigf, self.fphi])
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        self._args += ["-colin-fc", self.fphi.label()]
        self._args += ["-mtzout", "hklout.mtz"]

    def _result(self) -> ComitResult:
        self._check_files_exist("hklout.mtz")
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        return ComitResult(
            abcd=DataItem(mtz, "omit.ABCD"),
            fphi=DataItem(mtz, "omit.F_phi"),
            seconds=self._seconds,
        )
