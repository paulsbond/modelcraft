from typing import Optional, Union
import gemmi
from ..contents import AsuContents, PolymerType
from ..reflections import FsigF, FreeRFlag, ABCD, PhiFom, FPhi, write_mtz
from ..structure import write_mmcif
from .job import Job


class Parrot(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: FsigF,
        freer: FreeRFlag,
        phases: Union[ABCD, PhiFom],
        fphi: Optional[FPhi] = None,
        structure: Optional[gemmi.Structure] = None,
    ):
        super().__init__()
        args = []

        seqin = self.path("seqin.seq")
        args += ["-seqin", seqin]
        contents.write_sequence_file(seqin, PolymerType.PROTEIN)

        hklin = self.path("hklin.mtz")
        write_mtz(hklin, [fsigf, freer, phases, fphi])
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", freer.label()]
        if isinstance(phases, ABCD):
            args += ["-colin-hl", phases.label()]
        else:
            args += ["-colin-phifom", phases.label()]
        if fphi is not None:
            args += ["-colin-fc", fphi.label()]

        if structure is not None:
            xyzin = self.path("xyzin.cif")
            args += ["-pdbin", xyzin]  # or pdbin-ha or pdbin-mr?
            write_mmcif(xyzin, structure)

        args += ["-cycles", "5"]
        args += ["-anisotropy-correction"]

        hklout = self.path("hklout.mtz")
        args += ["-mtzout", hklout]

        self.run("cparrot", args)

        mtz = gemmi.read_mtz_file(hklout)
        self.abcd = ABCD(mtz, "parrot.ABCD")
        self.fphi = FPhi(mtz, "parrot.F_phi")

        self.finish()
