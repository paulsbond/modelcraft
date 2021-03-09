from math import isclose
from modelcraft.jobs import Acedrg
from modelcraft.residues import weight
import os
import gemmi


def test_0pr_smiles():
    smiles = "Cc1c(c(c(cn1)COP(=O)(O)O)CN[C@@H](Cc2ccc(cc2)O)C(=O)O)O"
    acedrg = Acedrg(code="0PR", smiles=smiles)
    weight = sum(atom.el.weight for atom in acedrg.chemcomp.atoms)
    assert isclose(weight, 409.31, abs_tol=0.01)
    acedrg.remove_files()


def test_gly_cif():
    path = os.path.join(os.environ["CLIBD"], "monomers", "g", "GLY.cif")
    cif = gemmi.cif.read(path)
    acedrg = Acedrg("GLY", cif=cif)
    old_weight = weight("GLY")
    new_weight = sum(atom.el.weight for atom in acedrg.chemcomp.atoms)
    assert isclose(old_weight, new_weight, abs_tol=0.01)
    acedrg.remove_files()
