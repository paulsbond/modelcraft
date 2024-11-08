# Changelog

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org).

## [5.0.0] - Unreleased

### Added

- EM --half-maps, --single-map and --build-map arguments.

### Removed

- EM --map and --blur arguments.

## [4.0.2] - 2024-10-24

### Fixed

- Crash when no residues are built.

## [4.0.1] - 2024-05-14

### Fixed

- Monomer library groups when a residue name is not in the library

## [4.0.0] - 2024-03-18

### Removed

- Coot RSR morph

## [3.6.0] - 2024-03-13

### Changed

- X-ray pipeline to run Buccaneer and Nautilus separately and combine the results
- Reporting numbers of protein and nucleic acid residues

## [3.5.3] - 2024-03-04

### Fixed

- Dummy atoms are now oxygen atoms instead of sodium atoms

## [3.5.2] - 2024-02-16

### Changed

- Output MTZ file to be the same as HKLOUT from the Refmac job

## [3.5.1] - 2024-01-18

### Changed

- Water chain naming to try to be compatible with PDB format

## [3.5.0] - 2024-01-09

### Added

- Option --threads to allow Buccaneer to use multiprocessing

## [3.4.0] - 2023-09-21

### Added

- Option --mask to provide a custom mask in EM mode instead of using EMDA mapmask

## [3.3.0] - 2023-08-01

### Changed

- Refinement now uses refmacat instead of refmac5

## [3.2.0] - 2023-07-10

### Changed

- Refine input model with and without Sheetbend and choose the best result
- Tests of monomer library entries
- TOXD pipeline test

## [3.1.1] - 2023-07-04

### Changed

- Refmac result to include initial R-work, R-free and FSC values
- Refmac test

## [3.1.0] - 2023-03-29

### Added

- Command line option: --overwrite-directory
- Writing current.cif and current.mtz while running

## [3.0.0] - 2023-03-09

### Added

- New EM pipeline starting from either halfmaps or a single map with optional blurring
- EMDA map mask job
- Servalcat trim, refine, FSC and NE map jobs
- Refmac map to MTZ job
- Libg job (not yet exposed)
- Support for a user-supplied restraints dictionary for ligands
- CCP-EM program tests
- Renaming SUL residues to SO4 and HOH O1 atoms to O on reading a structure

### Changed

- Renamed RefmacXray to Refmac
- Flushing print statements so they are written when piping to a file

### Removed

- Refmac EM job

## [2.4.1] - 2022-07-12

### Fixed

- Catching connection errors when requesting PDB entry contents
- Ensuring MTZ data items use the same ASU definition

## [2.4.0] - 2022-06-01

### Changed

- Comparing per-cycle outputs using R-work instead of R-free

### Fixed

- Nautilus test that was looking for a space in a residue name

## [2.3.6] - 2022-05-30

### Fixed

- Trimming spaces from residue names when reading structures

## [2.3.5] - 2022-05-23

### Changed

- Aligned columns in mmCIF output
- CCP4 map setup in EM mode for Gemmi version 0.5.3+

## [2.3.4] - 2022-05-19

### Fixed

- Consecutive residues function to yield last section for known structure IDs

## [2.3.3] - 2022-05-11

### Added

- Version and arguments in the output JSON file

### Fixed

- No longer using $CLIBD_MON/windows_reserved_words.txt

## [2.3.2] - 2022-03-11

### Fixed

- Buccaneer known structure function for modifications at N-terminus or joined together

## [2.3.1] - 2022-03-10

### Changed

- Reverted the 2.2.1 change to expected R-factors and FSC in the Refmac test

### Fixed

- Updating the cell when it is not needed or a chain has fewer than 3 atoms
- Ligands with CA atoms are now part of the known structure in Buccaneer

## [2.3.0] - 2022-02-24

### Added

- Command line arguments to disable individual steps
- Functions and classes available at the top package level
- Job information to the log and JSON file
- Residue and fragment metrics from Nautilus

### Fixed

- Monomer library paths for windows reserved words

## [2.2.3] - 2022-01-20

### Fixed

- Executable paths for Windows batch files (e.g. acedrg and coot)

## [2.2.2] - 2022-01-20

### Fixed

- Removed the dependency on CCP4i2 demo data from Molrep tests

## [2.2.1] - 2022-01-19

### Fixed

- CCP4i2 demo data path for Molrep tests in CCP4 8.0
- Expected R-factors and FSC for CCP4 8.0 Refmac test

## [2.2.0] - 2022-01-10

### Added

- Function to calculate geometry RMSZ for residues in the monomer library
- Argument tests

### Changed

- Runs inside a new directory to remove the chance of overwriting files
- Read multiple observation types (e.g. F,SIGF and I,SIGI) at once by default
- Errors if 0% or >50% of the reflections are in the free set

### Fixed

- Input phase labels for parrot and comit
- Main module execution using -m

## [2.1.0] - 2022-01-03

### Added

- Phase match job and tests
- Refmac resolution and completeness

### Fixed

- Termination after Buccaneer builds no residues

## [2.0.5] - 2021-12-03

### Added

- Usage of alpha and gamma angles when determining the spacegroup of a structure

## [2.0.4] - 2021-12-03

### Fixed

- Check on identical space groups due to differences such as R 3 vs H 3

## [2.0.3] - 2021-12-02

### Fixed

- Test that expected 0PR not to be in the monomer library.

## [2.0.2] - 2021-11-16

### Fixed

- Reading structures with gemmi versions between 0.4.8 and 0.5.0.

## [2.0.1] - 2021-11-15

### Added

- Cell module test.

### Fixed

- Removing atoms not in the monomer library after Nautilus.
- Updated PIP requirements.

## [2.0.0] - 2021-11-04

### Added

- Support for cryo-EM as well as X-ray data.
- Command line documentation.
- Printing of Refmac results to the log file.
- Generic checks that output files exist after running jobs.
- MOLREP job and tests.
- FREERFLAG job and test.
- Coot RSR morph job and test.
- Function to contract common column labels.
- Function to remove residues by name.
- Support for obsolete entries in modelcraft-contents.
- Decorator to run test functions in a temporary directory.
- Function to download files from PDBe
- Nautilus test
- Cryo-EM test
- Refmac jelly-body restraints option
- Sheetbend regularise option

### Changed

- Command line arguments.
- Ignoring file extensions when reading structures.
- Input phases not overwritten by dummy atom phases in first cycle.
- Cycles into a list in the JSON output.
- Setting model cell to the same as the MTZ file.

### Removed

- Polymer start parameter in the ASU contents description.
- Support for custom Buccaneer, Parrot and Sheetbend paths.

### Fixed

- Bug in the Coot side chains script.

## [1.0.0] - 2021-07-08

First non-pre-release.
