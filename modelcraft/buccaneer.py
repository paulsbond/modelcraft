

class Buccaneer(Job):
    def __init__(self, cycle, mtzin, pdbin=None, cycles=2):
        super().__init__(cycle, "Buccaneer")
        self.xmlout = self.path("xmlout.xml")
        self.xyzout = self.path("xyzout.pdb")
        self.run(mtzin, pdbin, cycles)
        self.finish()

    def run(self, mtzin, pdbin, cycles):
        stdin = []
        stdin.append("seqin " + args.seqin)
        stdin.append("colin-fo " + args.colin_fo)
        stdin.append("colin-free " + args.colin_free)
        if mtzin == args.mtzin:
            stdin.append("mtzin " + mtzin)
            if args.colin_hl is not None:
                stdin.append("colin-hl " + args.colin_hl)
            if args.colin_phifom is not None:
                stdin.append("colin-phifom " + args.colin_phifom)
        else:
            execute("cmtzjoin", [
                "-mtzout", self.path("hklin.mtz"),
                "-mtzin", args.mtzin, "-colin", args.colin_fo, "-colout", args.colin_fo,
                "-mtzin", args.mtzin, "-colin", args.colin_free, "-colout", args.colin_free,
                "-mtzin", mtzin, "-colin", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
                                 "-colout", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
                "-mtzin", mtzin, "-colin", "FWT,PHWT", "-colout", "FWT,PHWT",
            ])
            stdin.append("mtzin " + self.path("hklin.mtz"))
            stdin.append("colin-hl HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB")
            stdin.append("colin-fc FWT,PHWT")
        if pdbin is not None:
            stdin.append("pdbin " + pdbin)
            for structure in args.known_structure:
                stdin.append("known-structure " + structure)
        if args.mr_mode > 1:
            stdin.append("pdbin-mr " + args.mr_model)
            if args.mr_mode == 3:
                stdin.append("mr-model")
            if args.mr_mode in (4, 6):
                stdin.append("mr-model-filter")
            if args.mr_mode in (5, 6):
                stdin.append("mr-model-seed")
        stdin.append("pdbout " + self.xyzout)
        stdin.append("xmlout " + self.xmlout)
        stdin.append("anisotropy-correction")
        stdin.append("fast")
        stdin.append("cycles %d" % cycles)
        stdin.append("resolution 2.0")
        stdin.append("model-filter")
        stdin.append("model-filter-sigma 1.0")
        if args.semet:
            stdin.append("build-semet")
        execute(args.buccaneer, ["-stdin"], stdin)
        results = ET.parse(self.path("xmlout.xml")).getroot().find("Final")
        self.completeness_by_residues = float(results.find("CompletenessByResiduesBuilt").text)
        self.completeness_by_chains = float(results.find("CompletenessByChainsBuilt").text)
        self.chains_built = int(results.find("ChainsBuilt").text)
        self.fragments_built = int(results.find("FragmentsBuilt").text)
        self.residues_unique = int(results.find("ResiduesUnique").text)
        self.residues_built = int(results.find("ResiduesBuilt").text)
        self.residues_sequenced = int(results.find("ResiduesSequenced").text)
        self.longest_fragment = int(results.find("ResiduesLongestFragment").text)
