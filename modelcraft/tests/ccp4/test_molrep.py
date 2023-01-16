import os
import gemmi
from modelcraft.jobs.molrep import Molrep
from modelcraft.reflections import DataItem
from modelcraft.structure import read_structure
from . import ccp4_path


def test_rnase():
    pdb_path = ccp4_path("examples", "rnase", "mr_mod.pdb")
    mtz_path = ccp4_path("examples", "rnase", "rnase25.mtz")
    structure = read_structure(pdb_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    fmean = DataItem(mtz, "FNAT,SIGFNAT")
    molrep = Molrep(observations=fmean, structure=structure, number_of_monomers=2).run()
    assert molrep.err_level == 0
    assert molrep.mr_score > 0.3
    assert molrep.mr_zscore > 2.5
    assert molrep.n_solution == 2
    assert len(molrep.structure[0]) == 2


def test_hypf():
    directory = ccp4_path("examples", "mr_tutorial_2006", "data", "hypF")
    pdb_path = os.path.join(directory, "1v3z_B.pdb")
    mtz_path = os.path.join(directory, "hypF-1gxu-1gxt-HG_scaleit1.mtz")
    structure = read_structure(pdb_path)
    mtz = gemmi.read_mtz_file(mtz_path)
    fanom = DataItem(mtz, "F_Hg(+),SIGF_Hg(+),F_Hg(-),SIGF_Hg(-)")
    molrep = Molrep(observations=fanom, structure=structure).run()
    assert molrep.err_level == 0
    assert molrep.mr_score > 0.3
    assert molrep.mr_zscore > 2.5
    assert molrep.n_solution == 1
    assert len(molrep.structure[0]) == 1
