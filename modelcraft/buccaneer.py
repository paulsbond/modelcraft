from modelcraft.utils import run
from modelcraft.job import Job
import xml.etree.ElementTree as ET


class Buccaneer(Job):
    def __init__(self, args, directory, hklin, xyzin=None, cycles=2):
        super().__init__(directory)
        self.xmlout = self.path("xmlout.xml")
        stdin = self._get_stdin(args, hklin, xyzin, cycles)
        run(args.buccaneer, ["-stdin"], stdin, self.stdout, self.stderr)
        self._set_results()

    def _get_stdin(self, args, hklin, xyzin, cycles):
        stdin = []
        stdin.append("seqin %s" % args.seqin)
        stdin.append("colin-fo %s" % hklin.fsigf)
        stdin.append("colin-free %s" % args.colin_free)
        stdin.append("mtzin %s" % hklin.path)
        stdin.extend(self._colin_keywords(hklin))
        if xyzin is not None:
            stdin.append("pdbin %s" % xyzin)
            for structure in args.known_structure:
                stdin.append("known-structure %s" % structure)
        stdin.extend(self._mr_keywords(args))
        stdin.append("pdbout %s" % self.xyzout)
        stdin.append("xmlout %s" % self.xmlout)
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
        if args.mr_mode > 1:
            yield "pdbin-mr %s" % args.mr_model
            if args.mr_mode == 3:
                yield "mr-model"
            if args.mr_mode in (4, 6):
                yield "mr-model-filter"
            if args.mr_mode in (5, 6):
                yield "mr-model-seed"

    def _set_results(self):
        xml = ET.parse(self.xmlout).getroot().find("Final")
        self.completeness_by_residues = float(xml.find("CompletenessByResiduesBuilt").text)
        self.completeness_by_chains = float(xml.find("CompletenessByChainsBuilt").text)
        self.chains_built = int(xml.find("ChainsBuilt").text)
        self.fragments_built = int(xml.find("FragmentsBuilt").text)
        self.residues_unique = int(xml.find("ResiduesUnique").text)
        self.residues_built = int(xml.find("ResiduesBuilt").text)
        self.residues_sequenced = int(xml.find("ResiduesSequenced").text)
        self.longest_fragment = int(xml.find("ResiduesLongestFragment").text)
