class FindWaters(Job):
    def __init__(self, cycle, xyzin, hklin, dummy=False):
        super().__init__(cycle, "Finding waters")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, hklin, dummy)
        self.finish()

    def run(self, xyzin, hklin, dummy):
        arguments = [
            "--pdbin", xyzin,
            "--hklin", hklin,
            "--f", "FWT",
            "--phi", "PHWT",
            "--pdbout", self.path("waters.pdb"),
        ]
        if dummy:
            arguments.append("--flood")
        # --sigma 2.0
        # --min-dist X
        # --max-dist X
        # --flood-atom-radius 1.4 (adjusts contact distance)
        execute("findwaters", arguments)
        execute("pdb_merge", [
            "xyzin1", xyzin,
            "xyzin2", self.path("waters.pdb"),
            "xyzout", self.xyzout,
        ], ["NOMERGE", "END"])
