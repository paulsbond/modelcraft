import dataclasses
import gemmi
from ..job import Job
from ..pipeline import Pipeline
from ..reflections import DataItem


@dataclasses.dataclass
class RefmacResult:
    structure: gemmi.Structure
    abcd: DataItem
    fphi_best: DataItem
    fphi_diff: DataItem
    fphi_calc: DataItem
    rwork: float
    rfree: float


class Refmac(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem = None,
        cycles: int = 5,
        twinned: bool = False,
    ):
        super().__init__("refmac5")
        self._hklins["hklin.mtz"] = [fsigf, freer, phases]
        self._xyzins["xyzin.cif"] = structure
        self._args += ["HKLIN", "./hklin.mtz"]
        self._args += ["XYZIN", "./xyzin.cif"]
        self._args += ["HKLOUT", "./hklout.mtz"]
        self._args += ["XYZOUT", "./xyzout.cif"]
        self._args += ["XMLOUT", "./xmlout.xml"]
        labin = "FP=" + fsigf.label(0)
        labin += " SIGFP=" + fsigf.label(1)
        labin += " FREE=" + freer.label()
        if phases is not None:
            if phases.types == "AAAA":
                labin += " HLA=" + phases.label(0)
                labin += " HLB=" + phases.label(1)
                labin += " HLC=" + phases.label(2)
                labin += " HLD=" + phases.label(3)
            else:
                labin += " PHIB=" + phases.label(0)
                labin += " FOM=" + phases.label(1)
        self._stdin.append("LABIN " + labin)
        self._stdin.append("NCYCLES %d" % cycles)
        self._stdin.append("MAKE HYDR NO")
        if twinned:
            self._stdin.append("TWIN")
        self._stdin.append("MAKE NEWLIGAND NOEXIT")
        self._stdin.append("PHOUT")
        self._stdin.append("PNAME modelcraft")
        self._stdin.append("DNAME modelcraft")
        self._stdin.append("END")
        self._hklouts["hklout.mtz"] = None
        self._xmlouts["xmlout.xml"] = None
        self._xyzouts["xyzout.cif"] = None

    def run(self, pipeline: Pipeline = None) -> RefmacResult:
        super().run(pipeline)
        mtz = self._hklouts["hklout.mtz"]
        xml = self._xmlouts["xmlout.xml"]
        rworks = list(xml.iter("r_factor"))
        rfrees = list(xml.iter("r_free"))
        return RefmacResult(
            structure=self._xyzouts["xyzout.cif"],
            abcd=DataItem(mtz, "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fphi_calc=DataItem(mtz, "FC_ALL,PHIC_ALL"),
            rwork=float(rworks[-1].text) * 100,
            rfree=float(rfrees[-1].text) * 100,
        )
