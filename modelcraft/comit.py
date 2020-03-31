from modelcraft.reflections import DataFile
from modelcraft.job import Job


class Comit(Job):
    def __init__(self, directory, hklin):
        super().__init__(directory)
        self.run(
            "comit",
            ["-stdin"],
            [
                "mtzin " + hklin.path,
                "colin-fo " + hklin.fsigf,
                "colin-fc " + hklin.fcphic,
                "mtzout " + self.path("hklout.mtz"),
            ],
        )
        self.hklout = DataFile(self.path("hklout.mtz"))
        self.hklout.fsigf = hklin.fsigf
        self.hklout.free = hklin.free
        self.hklout.abcd = "omit.ABCD.A,omit.ABCD.B,omit.ABCD.C,omit.ABCD.D"
        self.hklout.fwphiw = "omit.F_phi.F,omit.F_phi.phi"
