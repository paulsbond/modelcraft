import gemmi
from modelcraft.data import DataItem, write_mtz
from modelcraft.job import Job
from modelcraft.model import write_mmcif


class Sheetbend(Job):
    def __init__(self, fsigf: DataItem, free: DataItem, structure: gemmi.Structure):
        super().__init__()
        args = []

        hklin = self.path("hklin.mtz")
        write_mtz(hklin, [fsigf, free])
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", free.label()]

        xyzin = self.path("xyzin.cif")
        write_mmcif(xyzin, structure)
        args += ["-pdbin", xyzin]

        xyzout = self.path("xyzout.cif")
        args += ["-pdbout", xyzout]

        args += ["-cycles", "12"]
        args += ["-resolution-by-cycle", "6,6,3"]

        self.run("csheetbend", args)
        self.structure = gemmi.read_structure(xyzout)
        self.finish()
