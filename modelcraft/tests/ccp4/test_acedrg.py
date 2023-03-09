import math
import os
import gemmi
from modelcraft.jobs.acedrg import Acedrg


def test_0pr_smiles():
    smiles = "Cc1c(c(c(cn1)COP(=O)(O)O)CN[C@@H](Cc2ccc(cc2)O)C(=O)O)O"
    acedrg = Acedrg(code="0PR", smiles=smiles).run()
    weight = sum(atom.el.weight for atom in acedrg.chemcomp.atoms)
    assert math.isclose(weight, 409.31, abs_tol=0.01)


def test_gly_cif():
    path = os.path.join(os.environ["CLIBD_MON"], "g", "GLY.cif")
    cif = gemmi.cif.read(path)
    acedrg = Acedrg("GLY", cif=cif).run()
    weight = sum(atom.el.weight for atom in acedrg.chemcomp.atoms)
    assert math.isclose(weight, 75.07, abs_tol=0.01)
