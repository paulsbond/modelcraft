residues_with_sidechains = {
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "MSE",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def neighbours_specs(residue):
    specs = [residue.spec]
    if residue.prev is not None:
        specs = [residue.prev.spec] + specs
    if residue.next is not None:
        specs = specs + [residue.next.spec]
    return specs


def refine(imol, specs):
    with AutoAccept():
        refine_residues_py(imol, specs)


def fix_side_chain(imol, imap, residue):
    spec = residue.spec
    delete_residue_sidechain(imol, spec[0], spec[1], spec[2], 0)
    neighbours = neighbours_specs(residue)
    refine(imol, neighbours)
    fill_partial_residue(imol, *spec)
    auto_fit_best_rotamer(spec[1], "", spec[2], spec[0], imol, imap, 1, 0.10)
    refine(imol, neighbours)


def fix_side_chains(imol, imap, imap_diff):
    model = Model(imol, imap, imap_diff)
    side_median = median([r.side_chain_correctness for r in model.residues if r.truncatable])
    threshold = side_median * 0.5
    for residue in model.residues:
        if (residue.name in residues_with_sidechains and not residue.truncatable) or (
            residue.truncatable and residue.side_chain_correctness < threshold
        ):
            fix_side_chain(imol, imap, residue)
