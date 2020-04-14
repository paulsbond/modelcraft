from modelcraft.coordinates import CoordinateFile
from modelcraft.job import Job


class Sheetbend(Job):
    def __init__(self, args, directory, xyzin):
        super().__init__(directory)
        arguments = []
        arguments += ["-mtzin", args.hklin.path]
        arguments += ["-colin-fo", args.hklin.fsigf]
        arguments += ["-colin-free", args.hklin.free]
        arguments += ["-pdbin", xyzin.path]
        arguments += ["-pdbout", self.path("xyzout.pdb")]
        arguments += ["-cycles", "12"]
        arguments += ["-resolution-by-cycle", "6,6,3"]
        self.run("csheetbend", arguments)
        self.xyzout = CoordinateFile(self.path("xyzout.pdb"))
        assert self.xyzout.residues == xyzin.residues
