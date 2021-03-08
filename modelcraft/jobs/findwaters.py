import gemmi
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif
from .job import Job


class FindWaters(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem, dummy: bool = False):
        name = "find_dummy" if dummy else "find_water"
        super().__init__(name)

        xyzin = self.path("xyzin.cif")
        hklin = self.path("hklin.mtz")
        water = self.path("water.pdb")  # In PDB format even with a cif extension

        write_mmcif(xyzin, structure)
        write_mtz(hklin, [fphi])

        args = []
        args += ["--pdbin", xyzin]
        args += ["--hklin", hklin]
        args += ["--f", fphi.label(0)]
        args += ["--phi", fphi.label(1)]
        if dummy:
            args += ["--flood"]
        args += ["--pdbout", water]
        self.run("findwaters", args)

        water_residues = []
        water_structure = read_structure(water)
        for water_chain in water_structure[0]:
            for water_residue in water_chain:
                water_residues.append(water_residue)

        self.structure = structure.clone()
        if len(water_residues) > 0:
            model = self.structure[0]
            chain = "DUM" if dummy else "WAT"
            if chain not in model:
                model.add_chain(chain)
            for residue in water_residues:
                model[chain].add_residue(residue)

        self.finish()
