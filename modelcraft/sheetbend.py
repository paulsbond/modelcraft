import gemmi
from modelcraft.data import DataItem, write_mtz
from modelcraft.job import Job
from modelcraft.model import write_mmcif


class Sheetbend(Job):
    def __init__(self, fsigf: DataItem, free: DataItem, structure: gemmi.Structure):
        super().__init__()

        hklin = self.path("hklin.mtz")
        xyzin = self.path("xyzin.cif")
        xyzout = self.path("xyzout.cif")

        write_mtz(hklin, [fsigf, free])
        write_mmcif(xyzin, structure)

        args = []
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", free.label()]
        args += ["-pdbin", xyzin]
        args += ["-pdbout", xyzout]
        args += ["-cycles", "12"]
        args += ["-resolution-by-cycle", "6,6,3"]
        self.run("csheetbend", args)

        self.structure = gemmi.read_structure(xyzout)

        self.finish()
