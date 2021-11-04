import os
import gemmi
from modelcraft.jobs.molrep import Molrep
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import ccp4_path


def test_gamma():
    directory = ccp4_path("share", "ccp4i2", "demo_data", "gamma")
    pdb_path = os.path.join(directory, "gamma_model.pdb")
    structure = read_structure(pdb_path)
    mtz_path = os.path.join(directory, "merged_intensities_Xe.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    ianom = DataItem(mtz, "Iplus,SIGIplus,Iminus,SIGIminus")
    molrep = Molrep(observations=ianom, structure=structure).run()
    assert molrep.err_level == 0
    assert molrep.n_solution == 1
    assert molrep.mr_score > 0.5
    assert molrep.mr_zscore > 5
    assert len(molrep.structure[0]) == molrep.n_solution


def test_gere():
    directory = ccp4_path("share", "ccp4i2", "demo_data", "gere")
    pdb_path = os.path.join(directory, "gere_molA.pdb")
    structure = read_structure(pdb_path)
    mtz_path = os.path.join(directory, "gere_scaled_data.mtz")
    mtz = gemmi.read_mtz_file(mtz_path)
    fmean = DataItem(mtz, "F_native,SIGF_native")
    molrep = Molrep(observations=fmean, structure=structure, number_of_monomers=6).run()
    assert molrep.err_level == 0
    assert molrep.n_solution >= 5
    assert molrep.mr_score > 0.5
    assert molrep.mr_zscore > 5
    assert len(molrep.structure[0]) == molrep.n_solution
