import dataclasses
from ..job import Job
from ..pipeline import Pipeline
from ..reflections import DataItem


@dataclasses.dataclass
class ComitResult:
    abcd: DataItem
    fphi: DataItem


class Comit(Job):
    def __init__(self, fsigf: DataItem, fphi: DataItem):
        super().__init__("comit")
        self._hklins["hklin.mtz"] = [fsigf, fphi]
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", fsigf.label()]
        self._args += ["-colin-fc", fphi.label()]
        self._args += ["-mtzout", "hklout.mtz"]
        self._hklouts["hklout.mtz"] = None

    def run(self, pipeline: Pipeline = None) -> ComitResult:
        super().run(pipeline)
        mtz = self._hklouts["hklout.mtz"]
        return ComitResult(
            abcd=DataItem(mtz, "omit.ABCD"),
            fphi=DataItem(mtz, "omit.F_phi"),
        )
