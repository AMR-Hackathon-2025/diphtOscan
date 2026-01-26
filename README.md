## _diphtOscan_

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

_diphtOscan_ is a command line script written in [Python](https://www.python.org/). _diphtOscan_ runs on UNIX, Linux and most OS X operating systems.
For more details, see the associated publication: [Hennart et al. 2023](https://www.biorxiv.org/content/10.1101/2023.02.20.529124v1).

_diphtOscan_ is a tool to search genomic assemblies of _Corynebacterium diphtheriae_ and other species of the _Corynebacterium diphtheriae_ species complex (CdSC) for:
* Species (e.g. _C. diphtheriae_, _C. ulcerans_, _C. pseudotuberculosis_, _C. belfantii_, _C. rouxii_ , _C. ramonii_ and _C. silvaticum_)
* Biovar-associated genes (_spuA_, nitrate reductase gene cluster)
* MLST sequence type
* Virulence factors, including _tox_ gene detection and disruption prediction
* Antimicrobial resistance determinants: acquired genes (_ermX_, _pbp2m_, …) and SNPs (e.g., _rpoB_, _gyrA_)
* Genomic context of genomic features associated with resistance
* Presence of integrons (using Integron Finder: https://github.com/gem-pasteur/Integron_Finder) 
* Tree building (using JolyTree: https://gitlab.pasteur.fr/GIPhy/JolyTree)

## Quick Start

```bash
# Install and setup
git clone https://github.com/AMR-Hackathon-2025/diphtOscan.git
cd diphtOscan
conda env create -f environment.yml
conda activate diphtoscan
pip install . --no-deps
diphtoscan update

# Run analysis (new subcommand format)
diphtoscan all -a your_assembly.fasta -st -t -res_vir -o results

# Or run individual analyses
diphtoscan species -a genome.fasta
diphtoscan mlst -a genome.fasta
diphtoscan amr -a genome.fasta
diphtoscan tox -a genome.fasta
diphtoscan qc -a genome.fasta
```

**Note:** The legacy flag-based syntax (e.g., `diphtoscan -a genome.fasta -st`) is still supported for backward compatibility, but the subcommand format is recommended.

## Installation and execution

To install:

1. Clone the [diphtoscan repository](https://github.com/AMR-Hackathon-2025/diphtOscan).
2. cd to the `diphtOscan` folder
3. Install the necessary dependencies using `conda` with `conda env create -f environment.yml`
4. Activate the `diphtoscan` environment with `conda activate diphtoscan`
5. Install the tool itself with `python -m pip install . --no-deps`
6. Update the database with `diphtoscan -u` before first using the tool.


## Usage

_diphtOscan_ uses a subcommand-based CLI structure. Run `diphtoscan --help` for an overview:

```
diphtoscan: a tool for characterising virulence and resistance in Corynebacterium

Available subcommands:
  all       Run complete analysis pipeline
  species   Species identification only
  mlst      MLST typing only
  amr       AMR and virulence gene screening
  tox       Tox allele detection only
  qc        Assembly QC metrics only
  update    Update MLST, Tox Allele & AMR databases
```

### Subcommand Examples

```bash
# Run complete analysis
diphtoscan all -a genome1.fasta genome2.fasta -o results

# Species identification only
diphtoscan species -a genome.fasta

# MLST typing only
diphtoscan mlst -a genome.fasta

# AMR/virulence screening
diphtoscan amr -a genome.fasta

# Tox allele detection
diphtoscan tox -a genome.fasta

# Assembly QC metrics only
diphtoscan qc -a genome.fasta

# Update databases
diphtoscan update
```

### Common Options

| Option | Description |
|--------|-------------|
| `-a, --assemblies` | FASTA file(s) for assemblies (required for analysis) |
| `-o, --outdir` | Output directory (default: `results_YYYY-MM-DD_HH-MM-SS_PID`) |
| `--format {tsv,json,html,all}` | Output format selection (default: `tsv`) |
| `-v, --verbose` | Enable verbose output with timestamps and debug information |
| `-q, --quiet` | Suppress all output except errors |
| `--no-qc` | Skip assembly quality control metrics calculation |
| `--min_identity` | Minimum alignment identity (default: 80) |
| `--min_coverage` | Minimum alignment coverage (default: 50) |
| `--threads` | Number of threads (default: 4) |
| `--overwrite` | Allow overwriting existing output directory |

### Analysis Options (for `all` subcommand)

| Option | Description |
|--------|-------------|
| `-st, --mlst` | Enable species and MLST sequence type detection |
| `-t, --tox` | Enable tox allele detection |
| `-res_vir, --resistance_virulence` | Enable resistance and virulence gene screening |
| `-plus, --extend_genotyping` | Enable extended virulence gene screening |
| `-integron` | Enable integron detection |
| `-tree` | Generate phylogenetic tree (requires >= 4 assemblies) |

### Legacy Syntax (Backward Compatible)

The original flag-based syntax is still supported:
```bash
diphtoscan -a genome.fasta -st -t -res_vir -o results
diphtoscan -u  # Update databases
```

## Example

In order to illustrate the usefulness of _diphtOscan_ and to describe its output files, the following use case example describes its usage for inferring a phylogenetic tree of _Corynebacterium diphtheriae_ genomes derived from the analysis of [Hennart et al](https://peercommunityjournal.org/articles/10.24072/pcjournal.307/).

##### Running _diphtOscan_

The following command line allows the script `diphtOscan` to be launched with default options on 8 threads:
```bash
# Using subcommand format (recommended)
diphtoscan all -a $genomes -st -res_vir --threads 8 -o Cdiphteriae

# With JSON output
diphtoscan all -a $genomes -st -res_vir --format all -o Cdiphteriae

# With HTML report only
diphtoscan all -a $genomes -st -res_vir --format html -o Cdiphteriae

# Verbose mode with detailed logging
diphtoscan all -a $genomes -st -res_vir -v -o Cdiphteriae

# Quiet mode (errors only)
diphtoscan all -a $genomes -st -res_vir -q -o Cdiphteriae
```

As the basename was set to 'Cdiphteriae', _diphtOscan_ writes the following output files:

* `Cdiphteriae.txt`: result file (tab-separated)
* `Cdiphteriae.json`: result file (JSON format, when `--format json` or `--format all`)
* `Cdiphteriae.html`: interactive HTML report (when `--format html` or `--format all`)
* `$strain.fa`: extracted sequences (for every assembly file)
* `$strain.out`: BLAST output file (for every assembly file)

**Note:** When processing multiple assemblies, a progress bar is displayed showing the processing status.

## Documentation

For detailed documentation, see the `docs/` folder:

- [Architecture Overview](docs/ARCHITECTURE.md) - System design and module responsibilities
- [Output Format](docs/OUTPUT_FORMAT.md) - Output files and annotation symbols explained
- [Algorithms](docs/ALGORITHMS.md) - Scientific methodology and thresholds
- [Database](docs/DATABASE.md) - Database structure and update procedures
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Common issues and FAQ

## Docker

diphtOscan is available as a Docker image for containerized execution.

### Building the Image

```bash
docker build -t diphtoscan .
```

### Running with Docker

```bash
# Show help
docker run --rm diphtoscan --help

# Run analysis (mount your data directory)
docker run --rm -v /path/to/data:/data -v /path/to/output:/output \
    diphtoscan all -a /data/assembly.fasta -st -t -res_vir -o /output/results

# Run with JSON output
docker run --rm -v /path/to/data:/data -v /path/to/output:/output \
    diphtoscan all -a /data/assembly.fasta -st --format json -o /output/results

# Update database inside container
docker run --rm -v diphtoscan-db:/opt/diphtoscan/data diphtoscan update
```

### Using Docker Compose

```yaml
version: '3'
services:
  diphtoscan:
    build: .
    volumes:
      - ./data:/data
      - ./output:/output
    command: ["all", "-a", "/data/assembly.fasta", "-st", "-o", "/output/results"]
```

## Running Tests

diphtOscan includes a comprehensive test suite with unit and integration tests.

```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=diphtoscan --cov-report=term-missing

# Run only unit tests
pytest tests/unit/ -v

# Run only integration tests
pytest tests/integration/ -v
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.

## Authors

- **Martin RETHORET-PASTY** - *Lead Developer* - Institut Pasteur
- Contact: martin.rethoret-pasty@pasteur.fr

For questions, bug reports, or feature requests, please open an issue on [GitHub](https://github.com/AMR-Hackathon-2025/diphtOscan/issues).

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

