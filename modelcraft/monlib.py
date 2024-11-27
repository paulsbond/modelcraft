import os
import gemmi
from .contents import PROTEIN_CODES, RNA_CODES, DNA_CODES


class MonLib(gemmi.MonLib):
    STANDARD: "MonLib"

    def __init__(self, resnames, libin: str = "", include_standard: bool = False):
        super().__init__()
        if libin:
            self.read_monomer_cif(libin)
        if include_standard:
            resnames = set(resnames)
            resnames |= set(PROTEIN_CODES.values())
            resnames |= set(RNA_CODES.values())
            resnames |= set(DNA_CODES.values())
            resnames |= {"HOH"}
        ok = self.read_monomer_lib(os.environ["CLIBD_MON"], list(resnames))
        if not ok:
            raise ValueError("Please create definitions for missing monomers.")

    def __contains__(self, code: str):
        return code in self.monomers

    def __getitem__(self, code: str):
        return self.monomers[code]

    def atom_ids(self, code: str):
        return {atom.id for atom in self[code].atoms}

    def group(self, code: str):
        return self[code].group if code in self else gemmi.ChemComp.Group.Null

    def is_nucleic(self, code: str) -> bool:
        return self.group(code) in {
            gemmi.ChemComp.Group.Dna,
            gemmi.ChemComp.Group.Rna,
            gemmi.ChemComp.Group.DnaRna,
        }

    def is_protein(self, code: str) -> bool:
        return self.group(code) in {
            gemmi.ChemComp.Group.Peptide,
            gemmi.ChemComp.Group.PPeptide,
            gemmi.ChemComp.Group.MPeptide,
        }

    def volume(self, code: str):
        return sum(18 for atom in self[code].atoms if not atom.is_hydrogen())

    def weight(self, code: str):
        return sum(atom.el.weight for atom in self[code].atoms)


MonLib.STANDARD = MonLib({}, include_standard=True)
