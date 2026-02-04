import dataclasses

import gemmi

from ..contents import AsuContents, PolymerType
from ..job import Job
from ..jobs.nautilus import NautilusResult
from ..reflections import DataItem, write_mtz
from ..structure import write_mmcif


@dataclasses.dataclass
class NucleoFindPrediction:
    phosphate: gemmi.Ccp4Map
    sugar: gemmi.Ccp4Map
    base: gemmi.Ccp4Map
    seconds: float


class NucleoFindPredict(Job):
    def __init__(self, fphi: DataItem):
        super().__init__("nucleofind")
        self.fphi = fphi

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fphi])
        self._args += ["--input", "hklin.mtz"]
        self._args += ["--amplitude", self.fphi.label(0)]
        self._args += ["--phase", self.fphi.label(1)]
        self._args += ["--model", "core"]
        self._args += ["--gpu"]
        self._args += ["--nthreads", "1"]
        self._args += ["--output", "."]

    def _result(self) -> NucleoFindPrediction:
        self._check_files_exist("nucleofind-phosphate.map")
        self._check_files_exist("nucleofind-sugar.map")
        self._check_files_exist("nucleofind-base.map")

        return NucleoFindPrediction(
            phosphate=gemmi.read_ccp4_map(self._path("nucleofind-phosphate.map")),
            sugar=gemmi.read_ccp4_map(self._path("nucleofind-sugar.map")),
            base=gemmi.read_ccp4_map(self._path("nucleofind-base.map")),
            seconds=self._seconds,
        )


class NucleoFindBuild(Job):
    def __init__(
        self,
        contents: AsuContents,
        fphi: DataItem,
        prediction: NucleoFindPrediction,
        structure: gemmi.Structure = None,
        em: bool = False,
    ):
        super().__init__("nucleofind-build")
        self.contents = contents
        self.fphi = fphi
        self.prediction = prediction
        self.structure = structure
        self.em = em

    def _setup(self) -> None:
        types = [PolymerType.RNA, PolymerType.DNA]
        self.contents.write_sequence_file(self._path("seqin.seq"), types)
        self._args += ["--seqin", "seqin.seq"]
        write_mtz(self._path("hklin.mtz"), [self.fphi])
        self._args += ["--mtzin", "hklin.mtz"]
        self._args += ["--colin-fc", self.fphi.label()]
        self.prediction.phosphate.write_ccp4_map(self._path("phosin.map"))
        self.prediction.sugar.write_ccp4_map(self._path("sugarin.map"))
        self.prediction.base.write_ccp4_map(self._path("basein.map"))
        self._args += ["--phosin", "phosin.map"]
        self._args += ["--sugarin", "sugarin.map"]
        self._args += ["--basein", "basein.map"]
        if self.structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.structure)
            self._args += ["--pdbin", "xyzin.cif"]
        self._args += ["--cycles", "3"]
        self._args += ["--pdbout", "xyzout.cif"]
        self._args += ["--xmlout", "xmlout.xml"]
        self._args += ["--remove_clashing_protein"]
        if self.em:
            self._args += ["--em"]

    def _result(self) -> NautilusResult:
        return NautilusResult(self)
