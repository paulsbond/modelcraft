from ...monlib import MonLib
from ...prune import _remove_isolated_residues, prune
from ...structure import ModelStats, read_structure
from . import ccp4_path, insulin_refmac


def test_insulin_prune():
    refmac = insulin_refmac()
    pruned = prune(
        structure=refmac.structure,
        fphi_best=refmac.fphi_best,
        fphi_diff=refmac.fphi_diff,
        fphi_calc=refmac.fphi_calc,
        residue_threshold=-2,
    )
    stats_in = ModelStats(refmac.structure)
    stats_out = ModelStats(pruned)
    assert stats_out.residues < stats_in.residues


def _chain_with_the_iron():
    "1rxf's randomised model: 182 waters, one iron and a 20-residue protein run."
    structure = read_structure(ccp4_path("examples", "data", "1rxf_randomise.pdb"))
    monlib = MonLib(structure[0].get_all_residue_names(), include_standard=True)
    return structure[0]["Y"], monlib


def _counts(chain, monlib):
    return (
        sum(1 for r in chain if r.name == "FE"),
        sum(1 for r in chain if r.name == "HOH"),
        sum(1 for r in chain if monlib.is_protein(r.name)),
    )


def test_ligands_and_waters_are_not_isolated_residues():
    chain, monlib = _chain_with_the_iron()
    assert _counts(chain, monlib) == (1, 182, 20)
    _remove_isolated_residues(chain, monlib)
    assert _counts(chain, monlib) == (1, 182, 20)


def test_an_isolated_residue_is_removed():
    chain, monlib = _chain_with_the_iron()
    protein = [i for i, r in enumerate(chain) if monlib.is_protein(r.name)]
    del chain[protein[11]]
    del chain[protein[9]]
    _remove_isolated_residues(chain, monlib)
    assert _counts(chain, monlib) == (1, 182, 17)
