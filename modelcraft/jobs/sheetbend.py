import dataclasses
import gemmi
from ..job import Job
from ..pipeline import Pipeline
from ..reflections import DataItem


@dataclasses.dataclass
class SheetbendResult:
    structure: gemmi.Structure


class Sheetbend(Job):
    def __init__(self, fsigf: DataItem, freer: DataItem, structure: gemmi.Structure):
        super().__init__("csheetbend")
        self._hklins["hklin.mtz"] = [fsigf, freer]
        self._xyzins["xyzin.cif"] = structure
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", fsigf.label()]
        self._args += ["-colin-free", freer.label()]
        self._args += ["-pdbin", "xyzin.cif"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cycles", "12"]
        self._args += ["-resolution-by-cycle", "6,3"]
        self._xyzouts["xyzout.cif"] = None

    def run(self, pipeline: Pipeline = None) -> SheetbendResult:
        super().run(pipeline)
        return SheetbendResult(structure=self._xyzouts["xyzout.cif"])
