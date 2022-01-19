import os
import gemmi
from modelcraft.jobs.molrep import Molrep
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import i2_demo_path


def test_gamma():
    pdb_path = i2_demo_path("gamma", "gamma_model.pdb")
    mtz_path = i2_demo_path("gamma", "merged_intensities_Xe.mtz")
    structure = read_structure(pdb_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    ianom = DataItem(mtz, "Iplus,SIGIplus,Iminus,SIGIminus")
    molrep = Molrep(observations=ianom, structure=structure).run()
    assert molrep.err_level == 0
    assert molrep.n_solution == 1
    assert molrep.mr_score > 0.5
    assert molrep.mr_zscore > 5
    assert len(molrep.structure[0]) == molrep.n_solution


def test_gere():
    pdb_path = i2_demo_path("gere", "gere_molA.pdb")
    mtz_path = i2_demo_path("gere", "gere_scaled_data.mtz")
    structure = read_structure(pdb_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    fmean = DataItem(mtz, "F_native,SIGF_native")
    molrep = Molrep(observations=fmean, structure=structure, number_of_monomers=6).run()
    assert molrep.err_level == 0
    assert molrep.n_solution >= 5
    assert molrep.mr_score > 0.5
    assert molrep.mr_zscore > 5
    assert len(molrep.structure[0]) == molrep.n_solution
