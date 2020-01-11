import Bio.SeqIO


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
            # TODO


class Polymer:
    def __init__(self, polymer_type, sequence, copies):
        self.polymer_type = polymer_type
        self.sequence = sequence
        self.copies = copies


class Protein(Polymer):
    codes = {
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "X",
        "Y",
        "Z",
    }

    def __init__(self, sequence, copies="Unknown"):
        self.polymer_type = "protein"


class Rna(Polymer):
    codes = {"A", "C", "G", "U", "X"}

    def __init__(self):
        self.polymer_type = "rna"


class Dna(Polymer):
    codes = {"A", "C", "G", "T", "X"}

    def __init__(self):
        self.polymer_type = "dna"


class Ligand:
    def __init__(self, code, copies="Unknown"):
        self.code = code
        self.copies = copies


class HeavyAtoms:
    def __init__(self, element, copies="Unknown"):
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
    unique_protein_codes = Protein.codes - Dna.codes
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
