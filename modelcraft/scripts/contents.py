import argparse
import sys
from ..contents import AsuContents, PolymerType


def main(argument_list=None):
    if argument_list is None:
        argument_list = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid")
    parser.add_argument("--contents-json", default="modelcraft-contents.json")
    parser.add_argument("--protein-fasta")
    parser.add_argument("--rna-fasta")
    parser.add_argument("--dna-fasta")
    args = parser.parse_args(argument_list)
    contents = AsuContents(args.pdbid)
    contents.write_json_file(args.contents_json)
    if args.protein_fasta is not None:
        contents.write_sequence_file(args.protein_fasta, PolymerType.PROTEIN)
    if args.rna_fasta is not None:
        contents.write_sequence_file(args.rna_fasta, PolymerType.RNA)
    if args.dna_fasta is not None:
        contents.write_sequence_file(args.dna_fasta, PolymerType.DNA)


if __name__ == "__main__":
    main()
