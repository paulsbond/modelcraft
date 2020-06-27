import gemmi
from .job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


class Sheetbend(Job):
    def __init__(self, fsigf: DataItem, freer: DataItem, structure: gemmi.Structure):
        super().__init__("sheetbend")
        args = []

        hklin = self.path("hklin.mtz")
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", freer.label()]
        write_mtz(hklin, [fsigf, freer])

        xyzin = self.path("xyzin.cif")
        args += ["-pdbin", xyzin]
        write_mmcif(xyzin, structure)

        args += ["-cycles", "12"]
        args += ["-resolution-by-cycle", "6,6,3"]

        xyzout = self.path("xyzout.cif")
        args += ["-pdbout", xyzout]

        self.run("csheetbend", args)

        self.structure = read_structure(xyzout)

        self.finish()
