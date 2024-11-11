import dataclasses
import gemmi
from .rscc import per_residue_rscc
from .contents import AsuContents, Polymer
from .reflections import DataItem


class _Fragment:
    def __init__(self, chain: gemmi.Chain, residue: gemmi.Residue, rsscs: dict):
        self.chain = chain
        self.residues = []
        self.rsccs = []
        self.fasta = ""
        self.protein = gemmi.find_tabulated_residue(residue.name).is_amino_acid()
        self.nucleic = gemmi.find_tabulated_residue(residue.name).is_nucleic_acid()
        self.add_residue(residue, rsscs)

    def addable(self, residue: gemmi.Residue) -> bool:
        info = gemmi.find_tabulated_residue(residue.name)
        previous = self.residues[-1]
        return (
            self.protein
            and info.is_amino_acid()
            and "C" in previous
            and "N" in residue
            and _min_dist(previous["C"], residue["N"]) < 2.0
        ) or (
            self.nucleic
            and info.is_nucleic_acid()
            and "O3'" in previous
            and "P" in residue
            and _min_dist(previous["O3'"], residue["P"]) < 2.0
        )

    def add_residue(self, residue: gemmi.Residue, rsscs: dict) -> None:
        self.residues.append(residue)
        self.rsccs.append(rsscs[(self.chain.name, str(residue.seqid))])
        if residue.name == "U" and "O4" not in residue:
            self.fasta += "X"
        else:
            self.fasta += gemmi.find_tabulated_residue(residue.name).fasta_code()

    def has_sequence(self):
        return self.fasta != "X" * len(self)

    def __len__(self):
        return len(self.residues)


def _fragments(model: gemmi.Model, rsccs: dict):
    for chain in model:
        fragment = _Fragment(chain, chain[0], rsccs)
        for residue in chain[1:]:
            if fragment.addable(residue):
                fragment.add_residue(residue, rsccs)
            else:
                yield fragment
                fragment = _Fragment(chain, residue, rsccs)
        if len(fragment.residues) > 0:
            yield fragment


def _min_dist(group1: gemmi.AtomGroup, group2: gemmi.AtomGroup) -> float:
    return min(atom1.pos.dist(atom2.pos) for atom1 in group1 for atom2 in group2)


class McModel:
    def __init__(
        self,
        structure: gemmi.Structure,
        contents: AsuContents,
        fphi: DataItem,
        model_index: int = 0,
    ):
        self.chains = []
        for polymer in contents.proteins + contents.rnas + contents.dnas:
            for _ in range(polymer.copies):
                self.chains.append(McChain(polymer))
        rsccs = per_residue_rscc(structure, fphi, radius=1.5, model_index=model_index)
        fragments = list(_fragments(structure[model_index], rsccs))
        for fragment in sorted(fragments, key=lambda f: sum(f.rsccs), reverse=True):
            if len(fragment) > 1 and fragment.has_sequence():
                print(
                    fragment.residues[0].seqid,
                    fragment.fasta,
                    fragment.residues[-1].seqid,
                )
                for chain in self.chains:
                    if chain.is_empty() or chain.name == fragment.chain.name:
                        chain.add(fragment)
                        break
                else:
                    overlaps = {}
                    for i, chain in enumerate(self.chains):
                        overlap = chain.overlap_info(fragment)
                        if overlap.linkable():
                            # trim using overlap
                            chain.add(fragment)
                            break
                        # Count contacts
        for chain in self.chains:
            print(chain.name)
            print(chain.residues)
        assert 1 == 0


@dataclasses.dataclass
class _OverlapInfo:
    slice_min: int
    slice_max: int
    linkable: bool
    contacts: int = 0


class McChain:
    def __init__(self, polymer: Polymer):
        self.name = None
        self.sequence = polymer.sequence
        self.residues = [None] * len(self)
        self.rsccs = [None] * len(self)

    def is_empty(self) -> bool:
        return self.residues == [None] * len(self)

    def add(self, fragment: _Fragment) -> None:
        if self.is_empty():
            self.name = fragment.chain.name
        for residue, rscc in zip(fragment.residues, fragment.rsccs):
            index = residue.seqid.num - 1
            self.residues[index] = residue.clone()
            self.rsccs[index] = rscc

    def overlap_info(self, fragment: _Fragment) -> _OverlapInfo:
        start_index = fragment.residues[0].seqid.num - 1
        slice_min = _slice_min(self.rsccs, fragment.rsccs, start_index)
        slice_max = len(fragment)  # TODO

    def __len__(self) -> float:
        return len(self.sequence)


def _slice_min(chain_rsccs: list, fragment_rsccs: list, start_index: int) -> int:
    len_ = len(fragment_rsccs)
    chain_rsccs = chain_rsccs[start_index : start_index + len_]
    chain_rsccs = [-1 if x is None else x for x in chain_rsccs]
    return max(range(len_ + 1), key=lambda i: sum(chain_rsccs[:i] + fragment_rsccs[i:]))
