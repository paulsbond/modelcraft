class Refmac(Job):
    def __init__(self, cycle, xyzin, use_phases):
        super().__init__(cycle, "REFMAC")
        self.hklout = self.path("hklout.mtz")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, use_phases)
        self.finish()

    def run(self, xyzin, use_phases):
        arguments = [
            "HKLIN", args.mtzin,
            "XYZIN", xyzin,
            "HKLOUT", self.hklout,
            "XYZOUT", self.xyzout,
            "XMLOUT", self.path("xmlout.xml"),
        ]
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
        stdin.append("NCYCLES 5")
        stdin.append("MAKE HYDR NO")
        if args.twinned:
            stdin.append("TWIN")
        stdin.append("MAKE NEWLIGAND NOEXIT")
        stdin.append("PHOUT")
        stdin.append("PNAME buccaneer")
        stdin.append("DNAME buccaneer")
        stdin.append("END")
        execute("refmac5", arguments, stdin)
        root = ET.parse(self.path("xmlout.xml")).getroot()
        rworks = list(root.iter("r_factor"))
        rfrees = list(root.iter("r_free"))
        self.initial_rwork = float(rworks[0].text)
        self.initial_rfree = float(rfrees[0].text)
        self.final_rwork = float(rworks[-1].text)
        self.final_rfree = float(rfrees[-1].text)
