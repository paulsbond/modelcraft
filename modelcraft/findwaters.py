from modelcraft.coordinates import CoordinateFile
from modelcraft.job import Job


class FindWaters(Job):
    def __init__(self, directory, xyzin, hklin, dummy=False):
        super().__init__(directory)
        arguments = [
            "--pdbin",
            xyzin.path,
            "--hklin",
            hklin.path,
            "--f",
            hklin.fphi.split(",")[0],
            "--phi",
            hklin.fphi.split(",")[0],
            "--pdbout",
            self.path("waters.pdb"),
        ]
        if dummy:
            arguments.append("--flood")
        # --sigma 2.0
        # --min-dist X
        # --max-dist X
        # --flood-atom-radius 1.4 (adjusts contact distance)
        self.run("findwaters", arguments)
        arguments = [
            "xyzin1",
            xyzin.path,
            "xyzin2",
            self.path("waters.pdb"),
            "xyzout",
            self.path("xyzout.pdb"),
        ]
        stdin = ["NOMERGE", "END"]
        self.run("pdb_merge", arguments, stdin)
        self.xyzout = CoordinateFile(self.path("xyzout.pdb"))
