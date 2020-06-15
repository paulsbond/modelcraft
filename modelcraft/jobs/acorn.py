# from .job import Job


# class Acorn(Job):
#     def __init__(self):
#         super().__init__()
#         self.finish()


# /home/paul/Programs/ccp4-7.0/bin/mtzdump \
# HKLIN "/home/paul/Desktop/tmp.mtz" \
# << EOF
# END
# EOF

# /home/paul/Programs/ccp4-7.0/bin/unique \
# HKLOUT "/tmp/paul/TEST_2_2_mtz.tmp" \
# << EOF
# cell 40.4450 42.6790 56.6510 90.0000 90.0000 90.0000
# symmetry 'P 21 21 21'
# resolution 1.0
# title [No title given]
# labout  F=F_unique SIGF=SIGF_unique
# end
# EOF

# /home/paul/Programs/ccp4-7.0/bin/cad \
# HKLIN1 "/home/paul/Desktop/tmp.mtz" \
# HKLIN2 "/tmp/paul/TEST_2_2_mtz.tmp" \
# HKLOUT "/tmp/paul/TEST_2_5_mtz.tmp" \
# << EOF
# LABIN FILE 1 ALL
# LABIN FILE 2 ALL
# END
# EOF

# /home/paul/Programs/ccp4-7.0/bin/mtzutils \
# HKLIN "/tmp/paul/TEST_2_5_mtz.tmp" \
# HKLOUT "/tmp/paul/TEST_2_7_mtz.tmp" \
# << EOF
# EXCL F_unique SIGF_unique
# END
# EOF

# /home/paul/Programs/ccp4-7.0/bin/phaser \
# << EOF
# MODE MR_ANO
# HKLIN "/tmp/paul/TEST_2_7_mtz.tmp"
# LABIN FP=FP SIGFP=SIGFP
# ROOT "/home/paul/Desktop/TEST_2"
# EOF

# /home/paul/Programs/ccp4-7.0/bin/ecalc \
# HKLIN "/home/paul/Desktop/TEST_2.mtz" \
# HKLOUT "/home/paul/Desktop/tmp_ecalc1.mtz" \
# << EOF
# title [No title given]
# labin  FP=FP_ISO SIGFP=SIGFP_ISO
# labout  FECALC=F_ECALC E=E_ISO SIGE=SIGE_ISO
# spacegroup 'P 21 21 21'
# EOF

# /home/paul/Programs/ccp4-7.0/bin/acorn \
# HKLIN "/home/paul/Desktop/tmp_ecalc1.mtz" \
# HKLOUT "/home/paul/Desktop/tmp_acorn1.mtz" \
# << EOF
# title [No title given]
# labin  FP=FP_ISO SIGFP=SIGFP_ISO E=E_ISO PHIN=hltofom.Phi_fom.phi WTIN=hltofom.Phi_fom.fom
# labout  PHIOUT=PHIOUT WTOUT=WTOUT

# #General ACORN keywords

# #Select Reflection Data
# EXTEND

# #Select Model Parameters

# #ACORN-MR keywords

# #ACORN-PHASE keywords
# NTRY 10

# END
# EOF

# /home/paul/Programs/ccp4-7.0/bin/fft \
# HKLIN "/home/paul/Desktop/tmp_acorn1.mtz" \
# MAPOUT "/tmp/paul/TEST_2_11_map.tmp" \
# << EOF
# title [No title given]
# labin -
#    F1=E_ISO SIG1=SIGE_ISO PHI=PHIOUT W=WTOUT
# end
# EOF

# /home/paul/Programs/ccp4-7.0/bin/mapmask \
# MAPIN "/tmp/paul/TEST_2_11_map.tmp" \
# MAPOUT "/tmp/paul/TEST_2.map" \
# << EOF
# SYMMETRY 'P 21 21 21'
# XYZLIM ASU
# EOF

# /home/paul/Programs/ccp4-7.0/bin/peakmax \
# MAPIN "/tmp/paul/TEST_2.map" \
# XYZOUT "/home/paul/Desktop/TEST_2.pdb" \
# XYZFRC "/home/paul/Desktop/TEST_2.ha" \
# << EOF
# threshhold -
#     rms -
#     5.0
# numpeaks 40
# output brookhaven frac
# EOF
