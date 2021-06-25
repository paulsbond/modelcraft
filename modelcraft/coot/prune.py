# CONSTS


protein_residues = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "MSE",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "UNK",
    "VAL",
}


main_chain_atoms = {" N  ", " CA ", " C  ", " O  ", " CB "}


# UTILS


def mean(values):
    return float(sum(values)) / len(values)


def median(values):
    count = len(values)
    if count < 1:
        return None
    i = count // 2
    if count % 2 == 1:
        return sorted(values)[i]
    return sum(sorted(values)[i - 1 : i + 1]) / 2.0


def distance(xyz1, xyz2):
    x_diff = xyz1[0] - xyz2[0]
    y_diff = xyz1[1] - xyz2[1]
    z_diff = xyz1[2] - xyz2[2]
    return (x_diff ** 2 + y_diff ** 2 + z_diff ** 2) ** 0.5


def is_protein(imol, name, res_spec):
    if name in protein_residues:
        return True
    atom_names = []
    atom_points = []
    for info in residue_info_py(imol, *res_spec):
        atom_names.append(info[0][0].strip())
        atom_points.append(info[2])
    if any(x not in atom_names for x in ("N", "CA", "C")):
        return False
    natoms = len(atom_names)
    n_xyz = [atom_points[i] for i in range(natoms) if atom_names[i] == "N"][0]
    ca_xyz = [atom_points[i] for i in range(natoms) if atom_names[i] == "CA"][0]
    c_xyz = [atom_points[i] for i in range(natoms) if atom_names[i] == "C"][0]
    return distance(n_xyz, ca_xyz) < 1.8 and distance(ca_xyz, c_xyz) < 1.8


# ATOM / RESIDUE / CHAIN / MODEL


class Atom:
    def __init__(self, atom_info):
        self.name = atom_info[0][0]
        self.point = atom_info[2]


class Residue:
    def __init__(self, model, spec, name):
        self.spec = spec
        self.name = name
        self.chain = spec[0]
        self.resno = spec[1]
        self.ins_code = spec[2]
        self.next = None
        self.prev = None
        self.atoms = {}
        self.truncatable = False
        for atom_info in residue_info_py(model.imol, *spec):
            atom = Atom(atom_info)
            model.atoms.append(atom)
            self.atoms[atom.name] = atom
            if atom.name not in main_chain_atoms:
                self.truncatable = True


class Chain:
    def __init__(self):
        self.residues = []
        self.correlations = []


class Model:
    def __init__(self, imol, imap):
        self.imol = imol
        self.imap = imap
        self.residues = []
        self.residue_spec_dict = {}
        self.atoms = []
        for ichain in range(n_chains(imol)):
            chain_id = chain_id_py(imol, ichain)
            for serial_num in range(chain_n_residues(chain_id, imol)):
                resno = seqnum_from_serial_number(imol, chain_id, serial_num)
                ins_code = insertion_code_from_serial_number(imol, chain_id, serial_num)
                spec = [chain_id, resno, ins_code]
                name = residue_name(imol, *spec)
                if is_protein(imol, name, spec):
                    residue = Residue(self, spec, name)
                    self.residues.append(residue)
                    self.residue_spec_dict[tuple(spec)] = residue
        self.set_connections()
        self.set_correlations()
        self.set_chains()

    def set_connections(self):
        for i in range(len(self.residues) - 1):
            res1 = self.residues[i]
            res2 = self.residues[i + 1]
            if " CA " not in res1.atoms or " C  " not in res1.atoms:
                continue
            if " N  " not in res2.atoms or " CA " not in res2.atoms:
                continue
            point1 = res1.atoms[" C  "].point
            point2 = res2.atoms[" N  "].point
            if distance(point1, point2) < 1.7:
                res1.next = res2
                res2.prev = res1

    def set_correlations(self):
        self.set_correlation("main_chain_correlation", 1)
        self.set_correlation("side_chain_correlation", 3)

    def set_correlation(self, attr, mask):
        specs = [residue.spec for residue in self.residues]
        for correlation in map_to_model_correlation_per_residue_py(
            self.imol, specs, mask, self.imap
        ):
            residue = self.residue_spec_dict[tuple(correlation[0][1:])]
            setattr(residue, attr, correlation[1])

    def set_chains(self):
        self.chains = {}
        for residue in self.residues:
            chain_id = residue.chain
            if chain_id not in self.chains:
                self.chains[chain_id] = Chain()
            chain = self.chains[chain_id]
            chain.residues.append(residue)
            if hasattr(residue, "main_chain_correlation"):
                chain.correlations.append(residue.main_chain_correlation)
        for chain in self.chains.values():
            if len(chain.correlations) > 0:
                chain.correlation = mean(chain.correlations)


# SCRIPTING


def prune(
    imol,
    imap,
    chains=True,
    chain_threshold="auto",
    max_chain_fraction=0.2,
    max_chain_length=20,
    residues=True,
    residue_threshold="auto",
    max_residue_fraction=0.2,
    remove_isolated_residues=True,
    sidechains=True,
    sidechain_threshold="auto",
    max_sidechain_fraction=0.2,
):

    model = Model(imol, imap)

    if len(model.chains) < 1:
        return

    if chains:
        main_median = median([r.main_chain_correlation for r in model.residues])
        if chain_threshold == "auto":
            chain_threshold = main_median * 0.2
        print(
            "PRUNE: Deleting chains (up to %d residues long) with scores < %.3f"
            % (max_chain_length, chain_threshold)
        )
        print(
            "PRUNE: Up to %.0f%% of residues will be deleted"
            % (max_chain_fraction * 100)
        )
        max_deleted = len(model.residues) * max_chain_fraction
        deleted = 0
        remaining = []
        chain_ids = model.chains.keys()
        chain_ids.sort(key=lambda cid: model.chains[cid].correlation)
        for chain_id in sorted(
            model.chains,
        ):
            chain = model.chains[chain_id]
            if (
                chain.correlation < chain_threshold
                and len(chain.residues) <= max_chain_length
                and deleted + len(chain.residues) <= max_deleted
            ):
                deleted += len(chain.residues)
                delete_chain(imol, chain_id)
            else:
                remaining.extend(chain.residues)
        print(
            "PRUNE: Deleted %.0f%% of residues"
            % (float(deleted) / len(model.residues) * 100)
        )
    else:
        remaining = model.residues

    if len(remaining) < 1:
        return

    if residues:
        main_median = median([r.main_chain_correlation for r in remaining])
        if residue_threshold == "auto":
            residue_threshold = main_median * 0.5
        print("PRUNE: Deleting residues with scores < %.3f" % residue_threshold)
        print(
            "PRUNE: Up to %.0f%% of residues will be deleted"
            % (max_residue_fraction * 100)
        )
        max_deleted = len(remaining) * max_residue_fraction
        deleted = 0
        remaining.sort(key=lambda r: r.main_chain_correlation)
        for residue in remaining:
            residue.delete = False
            if (
                residue.main_chain_correlation < residue_threshold
                and deleted + 1 <= max_deleted
            ):
                deleted += 1
                residue.delete = True
        if remove_isolated_residues:
            for residue in remaining:
                if residue.prev is not None and not residue.prev.delete:
                    continue
                if residue.next is not None and not residue.next.delete:
                    continue
                residue.delete = True
        for residue in remaining:
            if residue.delete:
                delete_residue(imol, residue.chain, residue.resno, residue.ins_code)
        print(
            "PRUNE: Deleted %.0f%% of residues"
            % (float(deleted) / len(remaining) * 100)
        )
        remaining = [r for r in remaining if not r.delete]

    remaining = [r for r in remaining if r.truncatable]
    if len(remaining) < 1:
        return

    if sidechains:
        side_median = median(
            [r.side_chain_correlation for r in model.residues if r.truncatable]
        )
        if sidechain_threshold == "auto":
            sidechain_threshold = side_median * 0.5
        print("PRUNE: Deleting sidechains with scores < %.3f" % sidechain_threshold)
        print(
            "PRUNE: Up to %.0f%% of sidechains will be deleted"
            % (max_sidechain_fraction * 100)
        )
        max_deleted = len(remaining) * max_sidechain_fraction
        deleted = 0
        remaining.sort(key=lambda r: r.side_chain_correlation)
        for residue in remaining:
            if (
                residue.side_chain_correlation < sidechain_threshold
                and deleted + 1 < max_deleted
            ):
                deleted += 1
                delete_residue_sidechain(
                    imol, residue.chain, residue.resno, residue.ins_code, 0
                )
        print(
            "PRUNE: Deleted %.0f%% of sidechains"
            % (float(deleted) / len(remaining) * 100)
        )
