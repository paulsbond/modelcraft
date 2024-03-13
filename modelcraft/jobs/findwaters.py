import dataclasses
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class FindWatersResult:
    structure: gemmi.Structure
    seconds: float


class FindWaters(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem, dummy: bool = False):
        super().__init__("findwaters")
        self.structure = structure
        self.fphi = fphi
        self.dummy = dummy

    def _setup(self) -> None:
        write_mtz(self._path("hklin.mtz"), [self.fphi])
        write_mmcif(self._path("xyzin.cif"), self.structure)
        self._args += ["--pdbin", "xyzin.cif"]
        self._args += ["--hklin", "hklin.mtz"]
        self._args += ["--f", self.fphi.label(0)]
        self._args += ["--phi", self.fphi.label(1)]
        if self.dummy:
            self._args += ["--flood"]
        self._args += ["--pdbout", "water.pdb"]

    def _result(self) -> FindWatersResult:
        self._check_files_exist("water.pdb")
        water_residues = []
        water_model = read_structure(self._path("water.pdb"))[0]
        for water_chain in water_model:
            for water_residue in water_chain:
                water_residue[0].element = gemmi.Element("O")
                water_residues.append(water_residue)
        structure = self.structure.clone()
        if len(water_residues) > 0:
            chain = gemmi.Chain("A")
            for residue in water_residues:
                chain.add_residue(residue)
            structure[0].add_chain(chain, unique_name=True)
        return FindWatersResult(
            structure=structure,
            seconds=self._seconds,
        )
