from modelcraft.coordinates import CoordinateFile
from modelcraft.reflections import ReflectionFile
from modelcraft.job import Job
import xml.etree.ElementTree as ET


class Refmac(Job):
    def __init__(self, args, directory, xyzin, cycles, use_phases=False):
        super().__init__(directory)
        self.hklout = ReflectionFile(
            path=self.path("hklout.mtz"),
            fsigf=args.hklin.fsigf,
            abcd="HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
            fphi="FWT,PHWT",
        )
        self.xmlout = self.path("xmlout.xml")
        arguments = self._get_arguments(args, xyzin)
        stdin = self._get_stdin(args, cycles, use_phases)
        self.run("refmac5", arguments, stdin)
        self.xyzout = CoordinateFile(self.path("xyzout.pdb"))
        self._set_results()

    def _get_arguments(self, args, xyzin):
        return [
            "HKLIN",
            args.hklin.path,
            "XYZIN",
            xyzin.path,
            "HKLOUT",
            self.hklout.path,
            "XYZOUT",
            self.path("xyzout.pdb"),
            "XMLOUT",
            self.xmlout,
        ]

    def _get_stdin(self, args, cycles, use_phases):
        stdin = []
        labin = "FP=" + args.colin_fp
        labin += " SIGFP=" + args.colin_sigfp
        labin += " FREE=" + args.colin_free
        if use_phases:
            if args.colin_hl is not None:
                labin += " HLA=" + args.colin_hla
                labin += " HLB=" + args.colin_hlb
                labin += " HLC=" + args.colin_hlc
                labin += " HLD=" + args.colin_hld
            if args.colin_phifom is not None:
                labin += " PHIB=" + args.colin_phi
                labin += " FOM=" + args.colin_fom
        stdin.append("LABIN " + labin)
        stdin.append("NCYCLES %d" % cycles)
        stdin.append("MAKE HYDR NO")
        if args.twinned:
            stdin.append("TWIN")
        stdin.append("MAKE NEWLIGAND NOEXIT")
        stdin.append("PHOUT")
        stdin.append("PNAME modelcraft")
        stdin.append("DNAME modelcraft")
        stdin.append("END")
        return stdin

    def _set_results(self):
        xml = ET.parse(self.xmlout).getroot()
        rworks = list(xml.iter("r_factor"))
        rfrees = list(xml.iter("r_free"))
        self.initial_rwork = float(rworks[0].text) * 100
        self.initial_rfree = float(rfrees[0].text) * 100
        self.final_rwork = float(rworks[-1].text) * 100
        self.final_rfree = float(rfrees[-1].text) * 100
