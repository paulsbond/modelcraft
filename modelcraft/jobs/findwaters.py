import gemmi
from ..reflections import DataItem, write_mtz
from ..structure import copy_structure, read_structure, write_mmcif
from .job import Job


class FindWaters(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem, dummy: bool = False):
        super().__init__()

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

        self.structure = copy_structure(structure)
        model = self.structure[0]
        chain = "dummy" if dummy else "water"
        if chain not in model:
            model.add_chain(chain)
        water_structure = read_structure(water)
        for water_chain in water_structure[0]:
            for water_residue in water_chain:
                model[chain].add_residue(water_residue)

        self.finish()
