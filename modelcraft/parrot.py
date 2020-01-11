from modelcraft.reflections import ReflectionFile
from modelcraft.job import Job


class Parrot(Job):
    def __init__(self, args, directory, hklin, xyzin=None):
        super().__init__(directory)
        hklin = self.create_hklin(args, hklin)
        self.hklout = ReflectionFile(
            path=self.path("hklout.mtz"),
            fsigf=hklin.fsigf,
            free=hklin.free,
            abcd="parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
            fphi="parrot.F_phi.F,parrot.F_phi.phi",
        )
        stdin = self._get_stdin(args, hklin, xyzin)
        self.run("cparrot", ["-stdin"], stdin)

    def _get_stdin(self, args, hklin, xyzin):
        stdin = []
        stdin.append("mtzin %s" % hklin.path)
        stdin.extend(self._colin_keywords(hklin))
        stdin.append("seqin %s" % args.seqin)
        if xyzin is not None:
            stdin.append("pdbin %s" % xyzin.path)  # or pdbin-ha or pdbin-mr?
        stdin.append("mtzout %s" % self.hklout.path)
        stdin.append("cycles 5")
        stdin.append("anisotropy-correction")
        return stdin

    def _colin_keywords(self, hklin):
        yield "colin-fo %s" % hklin.fsigf
        yield "colin-free %s" % hklin.free
        if hklin.abcd is not None:
            yield "colin-hl %s" % hklin.abcd
        if hklin.phifom is not None:
            yield "colin-phifom %s" % hklin.phifom
        if hklin.fphi is not None:
            yield "colin-fc %s" % hklin.fphi
