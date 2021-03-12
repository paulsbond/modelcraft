from typing import Optional
import gemmi
from ..contents import AsuContents
from ..job import Job
from ..polymer import PolymerType
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


class Nautilus(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: DataItem,
        freer: DataItem,
        phases: DataItem,
        fphi: Optional[DataItem] = None,
        structure: Optional[gemmi.Structure] = None,
        cycles: int = 3,
    ):
        super().__init__("nautilus")
        args = []

        seqin = self.path("seqin.seq")
        args += ["-seqin", seqin]
        contents.write_sequence_file(seqin, [PolymerType.RNA, PolymerType.DNA])

        hklin = self.path("hklin.mtz")
        data_items = [fsigf, freer, phases]
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", freer.label()]
        if phases.types == "AAAA":
            args += ["-colin-hl", phases.label()]
        else:
            args += ["-colin-phifom", phases.label()]
        if fphi is not None:
            args += ["-colin-fc", fphi.label()]
            data_items.append(fphi)
        write_mtz(hklin, data_items)

        if structure is not None:
            xyzin = self.path("xyzin.cif")
            args += ["-pdbin", xyzin]
            write_mmcif(xyzin, structure)

        args += ["-cycles", str(cycles)]
        args += ["-anisotropy-correction"]

        xyzout = self.path("xyzout.cif")
        args += ["-pdbout", xyzout]
        args += ["-cif"]

        self.run("cnautilus", args)
        self.structure = read_structure(xyzout)
        self.finish()
