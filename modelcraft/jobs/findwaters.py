import dataclasses
import gemmi
from ..job import Job
from ..pipeline import Pipeline
from ..reflections import DataItem


@dataclasses.dataclass
class FindWatersResult:
    structure: gemmi.Structure


class FindWaters(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem, dummy: bool = False):
        super().__init__("findwaters")
        self._structure = structure
        self._fphi = fphi
        self._dummy = dummy
        self._hklins["hklin.mtz"] = [fphi]
        self._xyzins["xyzin.cif"] = structure
        self._args += ["--pdbin", "xyzin.cif"]
        self._args += ["--hklin", "hklin.mtz"]
        self._args += ["--f", fphi.label(0)]
        self._args += ["--phi", fphi.label(1)]
        if dummy:
            self._args += ["--flood"]
        self._args += ["--pdbout", "water.pdb"]
        self._xyzouts["water.pdb"] = None

    def run(self, pipeline: Pipeline = None) -> FindWatersResult:
        super().run(pipeline)
        water_residues = []
        water_model = self._xyzouts["water.pdb"][0]
        for water_chain in water_model:
            for water_residue in water_chain:
                water_residues.append(water_residue)
        structure = self._structure.clone()
        if len(water_residues) > 0:
            model = structure[0]
            chain = "DUM" if self._dummy else "WAT"
            if chain not in model:
                model.add_chain(chain)
            for residue in water_residues:
                model[chain].add_residue(residue)
        return FindWatersResult(structure=structure)
