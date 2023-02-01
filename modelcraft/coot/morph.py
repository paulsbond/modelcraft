def rsr_morph(imol, imap, local_radius=5, gm_alpha=0.05, blur_b_factor=88):
    generate_self_restraints(imol, local_radius)
    set_show_extra_restraints(imol, 0)  # too confusing
    set_refinement_geman_mcclure_alpha(gm_alpha)
    set_draw_moving_atoms_restraints(1)  # not useful for non-graphics mode
    rmsd = map_sigma_py(imap)
    if rmsd is not None:
        imap_blurred = sharpen_blur_map(imap, blur_b_factor)
        set_imol_refinement_map(imap_blurred)
        set_matrix(15.0 / rmsd)
        residues = fit_protein_make_specs(imol, "all-chains")
        with AutoAccept():
            refine_residues_py(imol, residues)
