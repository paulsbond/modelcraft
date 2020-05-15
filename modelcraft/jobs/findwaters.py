import gemmi
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif
from .job import Job


class FindWaters(Job):
    def __init__(self, structure: gemmi.Structure, fphi: DataItem, dummy: bool = False):
        super().__init__()

        xyzin = self.path("xyzin.cif")
        hklin = self.path("hklin.mtz")
        waters = self.path("waters.pdb")  # In PDB format even with a cif extension
        xyzout = self.path("xyzout.cif")

        write_mmcif(xyzin, structure)
        write_mtz(hklin, [fphi])

        args = []
        args += ["--pdbin", xyzin]
        args += ["--hklin", hklin]
        args += ["--f", fphi.label(0)]
        args += ["--phi", fphi.label(1)]
        if dummy:
            # May be a bug in the minimum distance between the dummy atoms an structure
            args += ["--flood"]
        args += ["--pdbout", waters]
        self.run("findwaters", args)

        # TODO: Remove reliance on pdb_merge using gemmi
        args = []
        args += ["xyzin1", xyzin]
        args += ["xyzin2", waters]
        args += ["xyzout", xyzout]
        stdin = ["NOMERGE", "OUTPUT CIF", "END"]
        self.run("pdb_merge", args, stdin)

        self.structure = read_structure(xyzout)

        self.finish()
