import gemmi
from ..reflections import DataItem, write_mtz
from .job import Job


class Comit(Job):
    def __init__(self, fsigf: DataItem, fphi: DataItem):
        super().__init__("comit")
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
        self.abcd = DataItem(mtz, "omit.ABCD")
        self.fphi = DataItem(mtz, "omit.F_phi")

        self.finish()
