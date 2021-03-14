import dataclasses
import gemmi
from ..contents import AsuContents
from ..job import Job
from ..pipeline import Pipeline
from ..reflections import DataItem
from ..solvent import solvent_fraction


@dataclasses.dataclass
class ParrotResult:
    abcd: DataItem
    fphi: DataItem


class Parrot(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        structure: gemmi.Structure = None,
    ):
        super().__init__("cparrot")
        self._hklins["hklin.mtz"] = [fsigf, freer, phases, fphi]
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", fsigf.label()]
        self._args += ["-colin-free", freer.label()]
        if phases.types == "AAAA":
            self._args += ["-colin-hl", phases.label()]
        else:
            self._args += ["-colin-phifom", phases.label()]
        if fphi is not None:
            self._args += ["-colin-fc", fphi.label()]
        if structure is not None:
            self._xyzins["xyzin.cif"] = structure
            self._args += ["-pdbin-mr", "xyzin.cif"]
        self._args += ["-solvent-content", "%.3f" % solvent_fraction(contents, fsigf)]
        self._args += ["-cycles", "5"]
        self._args += ["-anisotropy-correction"]
        self._args += ["-mtzout", "hklout.mtz"]
        self._hklouts["hklout.mtz"] = None

    def run(self, pipeline: Pipeline = None) -> ParrotResult:
        super().run(pipeline)
        mtz = self._hklouts["hklout.mtz"]
        return ParrotResult(
            abcd=DataItem(mtz, "parrot.ABCD"),
            fphi=DataItem(mtz, "parrot.F_phi"),
        )
