import Bio.SeqIO
import modelcraft.residues as residues


class AsuContents:
    def __init__(self, relative=False):
        self.polymers = []
        self.ligands = []
        self.heavy_atoms = []
        self.relative = relative

    def add_polymers_from_sequence_file(self, path, polymer_type="auto"):
        for record in Bio.SeqIO.parse(path, "fasta"):
            sequence = str(record.seq).upper()
            if polymer_type == "auto":
                polymer_type = determine_polymer_type_from_sequence(sequence)
            if polymer_type == "protein":
                self.polymers.append(Protein(sequence))
            elif polymer_type == "rna":
                self.polymers.append(Rna(sequence))
            elif polymer_type == "dna":
                self.polymers.append(Dna(sequence))
            else:
                raise ValueError("Unknown polymer type: %s" % polymer_type)


class Polymer:
    def __init__(self, polymer_type, sequence, copies):
        self.polymer_type = polymer_type
        self.sequence = sequence
        self.copies = copies


class Protein(Polymer):
    def __init__(self, sequence, copies="unknown"):
        super().__init__("protein", sequence, copies)


class Rna(Polymer):
    def __init__(self, sequence, copies="unknown"):
        super().__init__("rna", sequence, copies)


class Dna(Polymer):
    def __init__(self, sequence, copies="unknown"):
        super().__init__("dna", sequence, copies)


class Ligand:
    def __init__(self, code, copies="unknown"):
        self.code = code
        self.copies = copies


class HeavyAtoms:
    def __init__(self, element, copies="unknown"):
        self.element = element
        self.copies = copies


def determine_polymer_type_from_sequence(sequence):
    def guess(polymer_type, specific=None):
        print("Guessing the polymer type of the following sequence:")
        print_sequence(sequence)
        print(
            "It is assumed to be %s" % (polymer_type if specific is None else specific)
        )
        return polymer_type

    codes = set(sequence)
    if "U" in codes:
        return "rna"
    protein_codes = {residue.code1 for residue in residues.protein}
    dna_codes = {residue.code1 for residue in residues.dna}
    unique_protein_codes = protein_codes - dna_codes
    if codes.intersection(unique_protein_codes):
        return "protein"
    if codes == {"A"}:
        return guess("protein", "polyalanine")
    if codes == {"G"}:
        return guess("protein", "polyglycine")
    if "T" in codes:
        return guess("dna")
    else:
        return guess("rna")


def print_sequence(sequence, line_length=60):
    for i in range(0, len(sequence), line_length):
        print(sequence[i : i + line_length])
