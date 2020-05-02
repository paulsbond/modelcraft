import gemmi
from ..reflections import FsigF, FPhi, ABCD, write_mtz
from .job import Job


class Comit(Job):
    def __init__(self, fsigf: FsigF, fphi: FPhi):
        super().__init__()
        args = []

        hklin = self.path("hklin.mtz")
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-fc", fphi.label()]
        write_mtz(hklin, [fsigf, fphi])

        hklout = self.path("hklout.mtz")
        args += ["-mtzout", hklout]

        self.run("comit", args)

        mtz = gemmi.read_mtz_file(hklout)
        self.abcd = ABCD(mtz, "omit.ABCD")
        self.fphi = FPhi(mtz, "omit.F_phi")

        self.finish()
