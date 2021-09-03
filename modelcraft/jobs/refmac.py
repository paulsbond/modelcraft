import dataclasses
import xml.etree.ElementTree as ET
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class RefmacResult:
    structure: gemmi.Structure
    abcd: DataItem
    fphi_best: DataItem
    fphi_diff: DataItem
    fphi_calc: DataItem
    rwork: float
    rfree: float


class _Refmac(Job):
    def __init__(self, structure: gemmi.Structure, cycles: int):
        super().__init__("refmac5")
        self.structure = structure
        self.cycles = cycles

    def _setup(self) -> None:
        write_mmcif(self._path("xyzin.cif"), self.structure)
        self._args += ["HKLIN", "./hklin.mtz"]
        self._args += ["XYZIN", "./xyzin.cif"]
        self._args += ["HKLOUT", "./hklout.mtz"]
        self._args += ["XYZOUT", "./xyzout.cif"]
        self._args += ["XMLOUT", "./xmlout.xml"]
        self._stdin.append("NCYCLES %d" % self.cycles)
        self._stdin.append("WEIGHT AUTO")
        self._stdin.append("MAKE HYDR NO")
        self._stdin.append("MAKE NEWLIGAND NOEXIT")
        self._stdin.append("PHOUT")
        self._stdin.append("PNAME modelcraft")
        self._stdin.append("DNAME modelcraft")
        self._stdin.append("END")

    def _result(self) -> RefmacResult:
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        xml = ET.parse(self._path("xmlout.xml")).getroot()
        rworks = list(xml.iter("r_factor"))
        rfrees = list(xml.iter("r_free"))
        return RefmacResult(
            structure=read_structure(self._path("xyzout.cif")),
            abcd=DataItem(mtz, "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fphi_calc=DataItem(mtz, "FC_ALL,PHIC_ALL"),
            rwork=float(rworks[-1].text) * 100,
            rfree=float(rfrees[-1].text) * 100,
        )


class RefmacXray(_Refmac):
    def __init__(
        self,
        structure: gemmi.Structure,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem = None,
        cycles: int = 5,
        twinned: bool = False,
    ):
        super().__init__(structure=structure, cycles=cycles)
        self.fsigf = fsigf
        self.freer = freer
        self.phases = phases
        self.twinned = twinned

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fsigf, self.freer, self.phases])
        labin = "FP=" + self.fsigf.label(0)
        labin += " SIGFP=" + self.fsigf.label(1)
        labin += " FREE=" + self.freer.label()
        if self.phases is not None:
            if self.phases.types == "AAAA":
                labin += " HLA=" + self.phases.label(0)
                labin += " HLB=" + self.phases.label(1)
                labin += " HLC=" + self.phases.label(2)
                labin += " HLD=" + self.phases.label(3)
            else:
                labin += " PHIB=" + self.phases.label(0)
                labin += " FOM=" + self.phases.label(1)
        self._stdin.append("LABIN " + labin)
        if self.twinned:
            self._stdin.append("TWIN")
        super()._setup()


class RefmacEm(_Refmac):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi: DataItem,
        freer: DataItem,
        cycles: int = 5,
    ):
        super().__init__(structure=structure, cycles=cycles)
        self.fphi = fphi
        self.freer = freer

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fphi, self.freer])
        labin = "FP=" + self.fphi.label(0)
        labin += " PHIB=" + self.fphi.label(1)
        labin += " FREE=" + self.freer.label()
        self._stdin.append("LABIN " + labin)
        self._stdin.append("SOURCE EM MB")
        self._stdin.append("SOLVENT NO")
        super()._setup()
