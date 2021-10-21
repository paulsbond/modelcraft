import dataclasses
import os
import gemmi
from ..contents import AsuContents, PolymerType
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class BuccaneerResult:
    structure: gemmi.Structure
    seconds: float


class Buccaneer(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        phases: DataItem,
        fphi: DataItem = None,
        freer: DataItem = None,
        input_structure: gemmi.Structure = None,
        mr_structure: gemmi.Structure = None,
        use_mr: bool = True,
        filter_mr: bool = True,
        seed_mr: bool = True,
        cycles: int = 2,
        em_mode: bool = False,
    ):
        super().__init__("cbuccaneer")
        self.contents = contents
        self.fsigf = fsigf
        self.phases = phases
        self.fphi = fphi
        self.freer = freer
        self.input_structure = input_structure
        self.mr_structure = mr_structure
        self.use_mr = use_mr
        self.filter_mr = filter_mr
        self.seed_mr = seed_mr
        self.cycles = cycles
        self.em_mode = em_mode

    def _setup(self) -> None:
        if self.em_mode:
            ref_dir = os.path.join(os.environ["CLIBD"], "reference_structures")
            ref_mtz = os.path.join(ref_dir, "reference-EMD-4116.mtz")
            ref_pdb = os.path.join(ref_dir, "reference-EMD-4116.pdb")
            self._args += ["-mtzin-ref", ref_mtz]
            self._args += ["-pdbin-ref", ref_pdb]
        types = [PolymerType.PROTEIN]
        self.contents.write_sequence_file(self._path("seqin.seq"), types)
        self._args += ["-seqin", "seqin.seq"]
        data_items = [self.fsigf, self.phases, self.fphi, self.freer]
        write_mtz(self._path("hklin.mtz"), data_items)
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", self.fsigf.label()]
        if self.phases.types == "AAAA":
            self._args += ["-colin-hl", self.phases.label()]
        else:
            self._args += ["-colin-phifom", self.phases.label()]
        if self.fphi is not None:
            self._args += ["-colin-fc", self.fphi.label()]
        if self.freer is not None:
            self._args += ["-colin-free", self.freer.label()]
        if self.input_structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.input_structure)
            self._args += ["-pdbin", "xyzin.cif"]
            self._args += ["-model-filter"]
            self._args += ["-model-filter-sigma", "1.0"]
            self._args += ["-nonprotein-radius", "1.0"]
        if self.mr_structure is not None:
            write_mmcif(self._path("xyzmr.cif"), self.mr_structure)
            self._args += ["-pdbin-mr", "xyzmr.cif"]
            if self.use_mr:
                self._args += ["-mr-model"]
                if self.filter_mr:
                    self._args += ["-mr-model-filter"]
                    self._args += ["-mr-model-filter-sigma", "2.0"]
                if self.seed_mr:
                    self._args += ["-mr-model-seed"]
        self._args += ["-cycles", str(self.cycles)]
        if self.contents.is_selenomet():
            self._args += ["-build-semet"]
        self._args += ["-fast"]
        self._args += ["-correlation-mode"]
        self._args += ["-anisotropy-correction"]
        self._args += ["-resolution", "2.0"]
        self._args += ["-pdbout", "xyzout.cif"]
        self._args += ["-cif"]

    def _result(self) -> BuccaneerResult:
        self._check_files_exist("xyzout.cif")
        return BuccaneerResult(
            structure=read_structure(self._path("xyzout.cif")),
            seconds=self._seconds,
        )
