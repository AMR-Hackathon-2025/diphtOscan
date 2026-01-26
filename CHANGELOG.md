# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Subcommand CLI structure**: New CLI with dedicated subcommands for focused analyses
  - `diphtoscan all` - Run complete analysis pipeline
  - `diphtoscan species` - Species identification only
  - `diphtoscan mlst` - MLST typing only
  - `diphtoscan amr` - AMR and virulence gene screening
  - `diphtoscan tox` - Tox allele detection only
  - `diphtoscan qc` - Assembly QC metrics only
  - `diphtoscan update` - Database update
- **JSON output format**: New `--format` option supporting `tsv`, `json`, `html`, or `all` output formats
  - Structured JSON with metadata, samples, and summary sections
  - Includes species counts and tox statistics in summary
- **HTML report generation**: Self-contained interactive HTML reports
  - Summary cards showing total samples, species distribution, tox prevalence, AMR prevalence
  - Searchable and sortable results table with color-coded cells
  - Species distribution chart with visual breakdown
  - Gene presence/absence heatmap coloring
  - No external dependencies - all CSS/JS embedded in single file
- **Verbose/quiet modes**: New `-v`/`--verbose` and `-q`/`--quiet` flags
  - Verbose mode (`-v`): Debug-level output with timestamps
  - Quiet mode (`-q`): Suppress all output except errors
  - Mutually exclusive options
- **Progress indicators**: tqdm-based progress bars for batch processing
  - Visual progress bar when processing multiple assemblies
  - Automatically disabled in quiet mode
- **Summary statistics**: End-of-analysis summary printed to console
  - Species distribution counts
  - Biovar distribution
  - Tox positive/negative counts with percentage
  - AMR prevalence by family
  - Assembly QC summary (average length, contig count)
- **Assembly QC module** (`assembly_qc.py`): Calculate quality metrics for input assemblies
  - Total genome length, contig count, N50/L50
  - GC content and ambiguous base (N) percentage
  - QC metrics included in results with `qc_` prefix
  - Skip with `--no-qc` flag
- **Comprehensive test suite**: 338 tests covering all modules
  - Unit tests for each module in `tests/unit/`
  - Integration tests for CLI in `tests/integration/`
  - Shared fixtures in `tests/conftest.py`
  - Test coverage reporting support
- **Docker support improvements**
  - Dockerfile for containerized execution
  - Docker Compose configuration example
  - Volume mounting for data and database persistence
- Comprehensive documentation in `docs/` directory
  - Architecture overview (ARCHITECTURE.md)
  - Output format documentation (OUTPUT_FORMAT.md)
  - Algorithm descriptions (ALGORITHMS.md)
  - Database documentation (DATABASE.md)
  - Troubleshooting guide (TROUBLESHOOTING.md)
- CONTRIBUTING.md with development guidelines and test documentation
- Quick Start section in README.md
- Badges for license and Python version in README.md
- NumPy-style docstrings to all major functions

### Changed
- CLI now uses subcommand structure (legacy flag syntax still supported)
- Updated README.md with subcommand examples and test instructions
- Updated README.md example to use `diphtoscan` command
- Fixed publication link in README.md

### Fixed
- Removed debug print statements from mlstBLAST.py

---

## [1.7.0] - 2024-08-21
### Changed
- Update tool `Integron_finder` to version `2.0.5`.
### Breaking Changes
- **Incompatibility with previous versions**: This update prevents the use of tool versions lower than `2.0.5`. Make sure you update the tool to version `2.0.5` or higher to continue using the project.

## [1.6.1] - 2024-08-21
### Added
- Completion of the Class & Subclass fields in the fam.tab file following the update of the AMRfinder database.
### Fixed
- Patch recurring error in dependency test function.
- Updated handling of deprecated pandas methods 

## [1.6.0] - 2024-08-20
### Fixed  
- Removal of sequences from the Corynebacterium_diphtheriae database following their addition to the AMRfinder database.
- Reallocation of "parent_node_id" after deletion of some due to AMRfinder database update.

## [1.5.0] - 2024-03-04
### Changed
- Formatted Coryne resistance database to AMRfinder 3.12 format.

---

[Unreleased]: https://github.com/AMR-Hackathon-2025/diphtOscan/compare/v1.7.0...HEAD
[1.7.0]: https://github.com/AMR-Hackathon-2025/diphtOscan/compare/v1.6.1...v1.7.0
[1.6.1]: https://github.com/AMR-Hackathon-2025/diphtOscan/compare/v1.6.0...v1.6.1
[1.6.0]: https://github.com/AMR-Hackathon-2025/diphtOscan/compare/v1.5.0...v1.6.0
[1.5.0]: https://github.com/AMR-Hackathon-2025/diphtOscan/releases/tag/v1.5.0