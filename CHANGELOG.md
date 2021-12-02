# Changelog

All notable changes to this project will be documented in this file.  
This project adheres to [Semantic Versioning](https://semver.org).

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
