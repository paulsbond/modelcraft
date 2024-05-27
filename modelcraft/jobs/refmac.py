import dataclasses
import shutil
import xml.etree.ElementTree as ET
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class RefmacResult:
    structure: gemmi.Structure
    mtz: gemmi.Mtz
    abcd: DataItem
    fphi_best: DataItem
    fphi_diff: DataItem
    fphi_calc: DataItem
    rwork: float
    rfree: float
    fsc: float
    initial_rwork: float
    initial_rfree: float
    initial_fsc: float
    data_completeness: float
    resolution_high: float
    seconds: float


class Refmac(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem = None,
        cycles: int = 5,
        twinned: bool = False,
        jelly_body: bool = False,
        libin: str = None,
    ):
        super().__init__("refmacat")
        self.structure = structure
        self.fsigf = fsigf
        self.freer = freer
        self.phases = phases
        self.cycles = cycles
        self.twinned = twinned
        self.jelly_body = jelly_body
        self.libin = libin

    def _setup(self) -> None:
        write_mmcif(self._path("xyzin.cif"), self.structure)
        write_mtz(self._path("hklin.mtz"), [self.fsigf, self.freer, self.phases])
        self._args += ["HKLIN", "./hklin.mtz"]
        self._args += ["XYZIN", "./xyzin.cif"]
        if self.libin:
            shutil.copy(self.libin, self._path("libin.cif"))
            self._args += ["LIBIN", "./libin.cif"]
        self._args += ["HKLOUT", "./hklout.mtz"]
        self._args += ["XYZOUT", "./xyzout.cif"]
        self._args += ["XMLOUT", "./xmlout.xml"]
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
        self._stdin.append(f"NCYCLES {self.cycles}")
        self._stdin.append("WEIGHT AUTO")
        if self.jelly_body:
            self._stdin.append("RIDGE DISTANCE SIGMA 0.02")
        if self.twinned:
            self._stdin.append("TWIN")
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
            mtz=mtz,
            abcd=DataItem(mtz, "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"),
            fphi_best=DataItem(mtz, "FWT,PHWT"),
            fphi_diff=DataItem(mtz, "DELFWT,PHDELWT"),
            fphi_calc=DataItem(mtz, "FC_ALL_LS,PHIC_ALL_LS"),
            rwork=float(rworks[-1].text),
            rfree=float(rfrees[-1].text),
            fsc=float(fscs[-1].text),
            initial_rwork=float(rworks[0].text),
            initial_rfree=float(rfrees[0].text),
            initial_fsc=float(fscs[0].text),
            data_completeness=float(xml.find("Overall_stats/data_completeness").text),
            resolution_high=float(xml.find("Overall_stats/resolution_high").text),
            seconds=self._seconds,
        )


@dataclasses.dataclass
class RefmacMapToMtzResult:
    fphi: DataItem
    seconds: float


class RefmacMapToMtz(Job):
    def __init__(self, density: gemmi.Ccp4Map, resolution: float, blur: float = 0):
        super().__init__("refmacat")
        self.density = density
        self.resolution = resolution
        self.blur = blur

    def _setup(self) -> None:
        self.density.write_ccp4_map(self._path("mapin.ccp4"))
        self._args += ["MAPIN", "mapin.ccp4"]
        self._args += ["HKLOUT", "hklout.mtz"]
        self._stdin.append("MODE SFCALC")
        self._stdin.append("SOURCE EM MB")
        self._stdin.append(f"RESOLUTION {self.resolution}")
        if self.blur > 0:
            self._stdin.append(f"SFCALC BLUR {self.blur}")
        elif self.blur < 0:
            self._stdin.append(f"SFCALC SHARP {-self.blur}")
        self._stdin.append("END")

    def _result(self) -> RefmacMapToMtzResult:
        self._check_files_exist("hklout.mtz")
        mtz = gemmi.read_mtz_file(self._path("hklout.mtz"))
        suffix = "0"
        if self.blur > 0:
            suffix = f"Blur_{self.blur:.2f}"
        elif self.blur < 0:
            suffix = f"Sharp_{-self.blur:.2f}"
        columns = f"Fout{suffix},Pout0"
        return RefmacMapToMtzResult(
            fphi=DataItem(mtz, columns),
            seconds=self._seconds,
        )
