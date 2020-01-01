from modelcraft.coordinates import CoordinateFile
from modelcraft.reflections import ReflectionFile
from modelcraft.job import Job


class Buccaneer(Job):
    def __init__(self, args, directory, hklin, xyzin=None, cycles=2, use_fphi=False):
        super().__init__(directory)
        hklin = self.create_hklin(args, hklin, use_fphi)
        stdin = self._get_stdin(args, hklin, xyzin, cycles)
        self.run(args.buccaneer, ["-stdin"], stdin)
        self.xyzout = CoordinateFile(self.path("xyzout.pdb"))

    def create_hklin(self, args, hklin, use_fphi):
        new_hklin = ReflectionFile(self.path("hklin.mtz"), "FP,SIGFP", "FREE", "HLA,HLB,HLC,HLD")
        args = [
            "-mtzout", self.path("hklin.mtz"),
            "-mtzin", args.hklin.path, "-colin", args.colin_fsigf, "-colout", "FP,SIGFP",
            "-mtzin", args.hklin.path, "-colin", args.colin_free, "-colout", "FREE",
            "-mtzin", hklin.path, "-colin", hklin.abcd, "-colout", "HLA,HLB,HLC,HLD",
        ]
        if use_fphi:
            args.extend(["-mtzin", hklin.path, "-colin", hklin.fphi, "-colout", "FWT,PHWT"])
            new_hklin.fphi = "FWT,PHWT"
        self.run("cmtzjoin", args)
        return new_hklin

    def _get_stdin(self, args, hklin, xyzin, cycles):
        stdin = []
        stdin.append("mtzin %s" % hklin.path)
        stdin.append("colin-fo %s" % hklin.fsigf)
        stdin.append("colin-free %s" % hklin.free)
        stdin.extend(self._colin_keywords(hklin))
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
        if hklin.abcd is not None:
            yield "colin-hl %s" % hklin.abcd
        if hklin.phifom is not None:
            yield "colin-phifom %s" % hklin.phifom
        if hklin.fphi is not None:
            yield "colin-fc %s" % hklin.fphi

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
