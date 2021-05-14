import dataclasses
import gemmi
from ..reflections import DataItem, write_mtz
from ..job import Job


@dataclasses.dataclass
class CTruncateResult:
    fmean: DataItem
    fanom: DataItem
    imean: DataItem


class CTruncate(Job):
    def __init__(self, observations: DataItem):
        super().__init__("ctruncate")
        self.observations = observations

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.observations])
        self._args += ["-hklin", "hklin.mtz"]
        self._args += ["-colin", f"/*/*/[{self.observations.label()}]"]
        self._args += ["-hklout", "hklout.mtz"]
        if self.observations.types == "KMKM":
            self._args += ["-Imean"]

    def _result(self) -> CTruncateResult:
        result = CTruncateResult(fmean=None, fanom=None, imean=None)
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        if self.observations.types == "KMKM":
            result.fmean = DataItem(mtz, "FMEAN,SIGFMEAN")
            result.fanom = DataItem(mtz, "F(+),SIGF(+),F(-),SIGF(-)")
            result.imean = DataItem(mtz, "_MEAN.I_sigI")
        elif self.observations.types == "JQ":
            result.fmean = DataItem(mtz, "F,SIGF")
        elif self.observations.types == "GLGL":
            result.fmean = DataItem(mtz, "FMEAN,SIGFMEAN")
        return result
