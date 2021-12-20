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
    fsc: float
    data_completeness: float
    resolution_high: float
    seconds: float


class _Refmac(Job):
    def __init__(self, structure: gemmi.Structure, cycles: int, jelly_body: bool):
        super().__init__("refmac5")
        self.structure = structure
        self.cycles = cycles
        self.jelly_body = jelly_body

    def _setup(self) -> None:
        write_mmcif(self._path("xyzin.cif"), self.structure)
        self._args += ["HKLIN", "./hklin.mtz"]
        self._args += ["XYZIN", "./xyzin.cif"]
        self._args += ["HKLOUT", "./hklout.mtz"]
        self._args += ["XYZOUT", "./xyzout.cif"]
        self._args += ["XMLOUT", "./xmlout.xml"]
        self._stdin.append("NCYCLES %d" % self.cycles)
        self._stdin.append("WEIGHT AUTO")
        if self.jelly_body:
            self._stdin.append("RIDGE DISTANCE SIGMA 0.02")
        self._stdin.append("MAKE HYDR NO")
        self._stdin.append("MAKE NEWLIGAND NOEXIT")
        self._stdin.append("PHOUT")
        self._stdin.append("PNAME modelcraft")
        self._stdin.append("DNAME modelcraft")
        self._stdin.append("END")

    def _result(self) -> RefmacResult:
        self._check_files_exist("xyzout.cif", "hklout.mtz", "xmlout.xml")
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        xml = ET.parse(self._path("xmlout.xml")).getroot()
        rworks = list(xml.iter("r_factor"))
        rfrees = list(xml.iter("r_free"))
        fscs = list(xml.iter("fscAver"))
        return RefmacResult(
            structure=read_structure(self._path("xyzout.cif")),
            abcd=DataItem(mtz, "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fphi_calc=DataItem(mtz, "FC_ALL,PHIC_ALL"),
            rwork=float(rworks[-1].text),
            rfree=float(rfrees[-1].text),
            fsc=float(fscs[-1].text),
            data_completeness=float(xml.find("Overall_stats/data_completeness").text),
            resolution_high=float(xml.find("Overall_stats/resolution_high").text),
            seconds=self._seconds,
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
        jelly_body: bool = False,
    ):
        super().__init__(structure=structure, cycles=cycles, jelly_body=jelly_body)
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
        cycles: int = 5,
        jelly_body: bool = False,
    ):
        super().__init__(structure=structure, cycles=cycles, jelly_body=jelly_body)
        self.fphi = fphi

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fphi])
        labin = "FP=" + self.fphi.label(0)
        labin += " PHIB=" + self.fphi.label(1)
        self._stdin.append("LABIN " + labin)
        self._stdin.append("SOURCE EM MB")
        self._stdin.append("SOLVENT NO")
        super()._setup()
