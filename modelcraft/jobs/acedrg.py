import dataclasses
import gemmi
from ..job import Job
from ..pipeline import Pipeline


@dataclasses.dataclass
class AcedrgResult:
    chemcomp: gemmi.ChemComp


class Acedrg(Job):
    def __init__(self, code: str, cif: gemmi.cif.Document = None, smiles: str = None):
        super().__init__("acedrg")
        self._args += ["-r", code.upper()]
        self._args += ["-o", "output"]
        if cif is not None:
            self._cifins["input.cif"] = cif
            self._args += ["-c", "input.cif"]
        if smiles is not None:
            self._args += ["-i", smiles]
        self._cifouts["output.cif"] = None

    def run(self, pipeline: Pipeline = None) -> AcedrgResult:
        super().run(pipeline)
        block = self._cifouts["output.cif"][-1]
        return AcedrgResult(chemcomp=gemmi.make_chemcomp_from_block(block))
