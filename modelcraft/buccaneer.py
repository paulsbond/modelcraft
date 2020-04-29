from typing import List, Optional, Union
import gemmi
from modelcraft.contents import AsuContents
from modelcraft.job import Job
from modelcraft.reflections import FsigF, FreeRFlag, ABCD, PhiFom, FPhi, write_mtz
from modelcraft.structure import write_mmcif


class Buccaneer(Job):
    def __init__(
        self,
        contents: AsuContents,
        fsigf: FsigF,
        freer: FreeRFlag,
        phases: Union[ABCD, PhiFom],
        fphi: Optional[FPhi] = None,
        input_structure: Optional[gemmi.Structure] = None,
        known_structure: Optional[List[str]] = None,
        mr_structure: Optional[gemmi.Structure] = None,
        use_mr: bool = True,
        filter_mr: bool = True,
        seed_mr: bool = True,
        cycles: int = 2,
        semet: bool = False,
        program: str = "cbuccaneer",
    ):
        super().__init__()
        args = []

        seqin = self.path("seqin.seq")
        args += ["-seqin", seqin]
        contents.write_protein(seqin)

        hklin = self.path("hklin.mtz")
        data_items = [fsigf, freer, phases]
        args += ["-mtzin", hklin]
        args += ["-colin-fo", fsigf.label()]
        args += ["-colin-free", freer.label()]
        if phases is ABCD:
            args += ["-colin-hl", phases.label()]
        else:
            args += ["-colin-phifom", phases.label()]
        if fphi is not None:
            args += ["-colin-fc", fphi.label()]
            data_items.append(fphi)
        write_mtz(hklin, data_items)

        if input_structure is not None:
            xyzin = self.path("xyzin.cif")
            args += ["-pdbin", xyzin]
            if known_structure is not None:
                for keyword in known_structure:
                    args += ["-known-structure", keyword]
            write_mmcif(xyzin, input_structure)

        if mr_structure is not None:
            xyzmr = self.path("xyzmr.cif")
            args += ["-pdbin-mr", xyzmr]
            if use_mr:
                args += ["mr-model"]
                if filter_mr:
                    args += ["mr-model-filter"]
                    args += ["mr-model-filter-sigma 2.0"]
                if seed_mr:
                    args += ["mr-model-seed"]
            write_mmcif(xyzmr, mr_structure)

        args += ["-cycles", cycles]
        if semet:
            args += ["-build-semet"]
        args += ["-fast"]
        args += ["-correlation-mode"]
        args += ["-anisotropy-correction"]
        args += ["-resolution 2.0"]
        args += ["-model-filter"]
        args += ["-model-filter-sigma 1.0"]

        xyzout = self.path("xyzout.cif")
        args += ["-pdbout", xyzout]

        self.run(program, args)
        self.structure = gemmi.read_structure(xyzout)
        self.finish()
