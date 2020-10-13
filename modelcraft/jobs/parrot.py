from typing import Optional
import gemmi
from ..contents import AsuContents, PolymerType
from ..reflections import DataItem, write_mtz
from ..structure import write_mmcif
from .job import Job


class Parrot(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem,
        fphi: Optional[DataItem] = None,
        structure: Optional[gemmi.Structure] = None,
    ):
        super().__init__("parrot")
        args = []

        seqin = self.path("seqin.seq")
        args += ["-seqin", seqin]
        contents.write_sequence_file(seqin, PolymerType.PROTEIN)

        hklin = self.path("hklin.mtz")
        write_mtz(hklin, [fsigf, freer, phases, fphi])
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", freer.label()]
        if phases.types == "AAAA":
            args += ["-colin-hl", phases.label()]
        else:
            args += ["-colin-phifom", phases.label()]
        if fphi is not None:
            args += ["-colin-fc", fphi.label()]

        if structure is not None:
            xyzin = self.path("xyzin.cif")
            args += ["-pdbin-mr", xyzin]
            write_mmcif(xyzin, structure)

        args += ["-cycles", "5"]
        args += ["-anisotropy-correction"]

        hklout = self.path("hklout.mtz")
        args += ["-mtzout", hklout]

        self.run("cparrot", args)

        mtz = gemmi.read_mtz_file(hklout)
        self.abcd = DataItem(mtz, "parrot.ABCD")
        self.fphi = DataItem(mtz, "parrot.F_phi")

        self.finish()
