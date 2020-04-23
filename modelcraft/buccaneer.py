import gemmi
from modelcraft.data import DataItem, write_mtz
from modelcraft.job import Job
from modelcraft.model import write_mmcif
from modelcraft.contents import AsuContents


class Buccaneer(Job):
    def __init__(
        self,
        fsigf: DataItem,
        free: DataItem,
        phases: DataItem,
        contents: AsuContents,
        fphi: DataItem = None,
        structure: gemmi.Structure = None,
        cycles: int = 2,
        program: str = "cbuccaneer",
    ):
        super().__init__()
        args = []

        hklin = self.path("hklin.mtz")
        data_items = [fsigf, free, phases]
        args += "mtzin %s" % hklin.path)
        write_mtz(hklin, data_items)
        args += ["colin-fo", fsigf.label()]
        args += ["colin-free", free.label()]
        if hklin.abcd is not None:
            yield "colin-hl %s" % hklin.abcd
        if hklin.phifom is not None:
            yield "colin-phifom %s" % hklin.phifom
        if hklin.fwphiw is not None:
            yield "colin-fc %s" % hklin.fwphiw

        xyzin = self.path("xyzin.cif")
        seqin = self.path("seqin.seq")
        xyzout = self.path("xyzout.cif")

        write_mmcif(xyzin, structure)

        self.run(program, args)
        self.structure = gemmi.read_structure(xyzout)
        self.finish()

        stdin.append("seqin %s" % args.seqin)
        if xyzin is not None:
            stdin.append("pdbin %s" % xyzin.path)
            for structure in args.known_structure:
                stdin.append("known-structure %s" % structure)
        stdin.extend(self._mr_keywords(args))
        stdin.append("pdbout %s" % self.path("xyzout.pdb"))
        stdin.append("cycles %d" % cycles)
        if args.semet:
            stdin.append("build-semet")
        stdin.append("fast")
        stdin.append("correlation-mode")
        stdin.append("anisotropy-correction")
        stdin.append("resolution 2.0")
        stdin.append("model-filter")
        stdin.append("model-filter-sigma 1.0")
        return stdin

    def _colin_keywords(self, hklin):


    def _mr_keywords(self, args):
        if args.mr_model is not None and args.mr_mode > 1:
            yield "pdbin-mr %s" % args.mr_model.path
            if args.mr_mode == 3:
                yield "mr-model"
            if args.mr_mode in (4, 6):
                yield "mr-model-filter"
                yield "mr-model-filter-sigma 2.0"
            if args.mr_mode in (5, 6):
                yield "mr-model-seed"
