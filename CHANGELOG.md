# Changelog

All notable changes to this project will be documented in this file.  
This project adheres to [Semantic Versioning](https://semver.org).

## [2.3.0] - Unreleased

### Added

- Command line arguments to disable individual steps
- Functions and classes available at the top package level
- Job information to the log and JSON file

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
