# diphtOscan Architecture

This document describes the system design and architecture of diphtOscan.

## Overview

diphtOscan is a bioinformatics tool for characterizing virulence and antimicrobial resistance in *Corynebacterium diphtheriae* and related species of the *C. diphtheriae* species complex (CdSC).

## Module Responsibilities

```
diphtoscan/
├── cli.py                  # Command-line interface and main orchestration
├── species.py              # Species identification using Mash distance
├── mlstBLAST.py            # MLST sequence type calling via BLAST
├── blastn.py               # BLAST wrapper and hit processing
├── truncation.py           # Gene truncation detection at amino acid level
├── assembly_qc.py          # Assembly quality control metrics
├── logging_config.py       # Logging configuration (verbose/quiet modes)
├── summary.py              # Summary statistics calculation and formatting
├── html_report.py          # HTML report generation using Jinja2
├── utils.py                # Utility functions, AMRFinder output processing
├── template_iTOL.py        # iTOL visualization template generation
├── updating_database.py    # Database update functionality
├── download_alleles_st.py  # MLST allele/profile downloading from PubMLST
├── jolytree_generation.py  # Phylogenetic tree generation using JolyTree
├── misc.py                 # Miscellaneous helper functions
└── templates/
    └── report.html         # Jinja2 template for HTML reports
```

### Module Descriptions

#### `cli.py`
The main entry point for diphtOscan. Responsibilities:
- Parse command-line arguments
- Validate dependencies (Mash, BLAST+, AMRFinderPlus, etc.)
- Orchestrate the analysis pipeline
- Manage output directory creation
- Coordinate results aggregation

#### `species.py`
Performs species identification using Mash genomic distance estimation.
- Computes Mash distances against reference species sketches
- Classifies matches as "strong" (<=0.05), "weak" (<=0.1), or "unknown"
- Determines if genome belongs to the *C. diphtheriae* species complex

#### `mlstBLAST.py`
Implements MLST (Multi-Locus Sequence Typing) calling.
- Runs BLAST searches against allele databases
- Matches alleles to sequence type (ST) profiles
- Handles imprecise matches and locus variants
- Reports novel alleles and partial matches

#### `blastn.py`
Provides BLAST nucleotide search functionality.
- Builds BLAST databases as needed
- Executes blastn searches with optimized parameters
- Parses BLAST output into structured `BlastHit` objects
- Culls redundant overlapping hits

#### `truncation.py`
Detects gene truncations at the amino acid level.
- Translates nucleotide sequences using bacterial genetic code
- Checks for premature stop codons
- Calculates amino acid coverage percentage

#### `assembly_qc.py`
Calculates assembly quality control metrics from input FASTA files.
- **Genome size**: Total assembly length in base pairs
- **Contig statistics**: Number of contigs, N50, L50, largest contig
- **GC content**: Percentage of G+C bases
- **N content**: Percentage of ambiguous bases (N)
- Provides `calculate_assembly_stats()` for comprehensive QC analysis
- Results are prefixed with `qc_` in output (e.g., `qc_n50`, `qc_gc_percent`)

#### `logging_config.py`
Provides centralized logging configuration for diphtOscan.
- Configures logging levels based on verbosity flags
- **Verbose mode** (`-v`): DEBUG level with timestamps
- **Quiet mode** (`-q`): ERROR level only
- **Default**: INFO level with minimal formatting
- Provides helper functions: `log_debug()`, `log_info()`, `log_warning()`, `log_error()`

#### `summary.py`
Generates summary statistics from batch analysis results.
- Calculates species distribution counts
- Computes biovar distribution
- Calculates tox prevalence (count and percentage)
- Summarizes AMR prevalence by gene family
- Aggregates QC metrics (average length, contigs)
- Formats summary for console output or JSON embedding

#### `html_report.py`
Generates self-contained interactive HTML reports using Jinja2 templates.
- Creates summary cards with visual statistics
- Produces searchable, sortable results tables
- Color-codes cells based on gene presence/absence
- Embeds all CSS and JavaScript for portability
- Uses template in `templates/report.html`

#### `utils.py`
Contains utility functions for data processing.
- Processes AMRFinderPlus output into tabular format
- Calculates genomic context for resistance genes
- Manages database paths and file operations
- Defines virulence gene categories

#### `template_iTOL.py`
Generates iTOL (Interactive Tree of Life) visualization files.
- Creates binary presence/absence datasets
- Generates color strip annotations
- Outputs toxin and AMR family visualizations

#### `updating_database.py`
Handles database updates from external sources.
- Downloads MLST profiles from PubMLST
- Fetches AMRFinderPlus database files
- Merges custom *C. diphtheriae* annotations

#### `download_alleles_st.py`
Downloads allele sequences and ST profiles from PubMLST API.
- Fetches allele FASTA files for each locus
- Downloads ST profile definitions
- Creates concatenated sequence databases

#### `jolytree_generation.py`
Generates phylogenetic trees using JolyTree.
- Prepares input assembly files
- Executes JolyTree pipeline
- Manages temporary working directories

## Data Flow

```
Input FASTA Assembly(s)
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│               Assembly QC Metrics (unless --no-qc)               │
│  (assembly_qc.py → genome size, N50, GC content, etc.)          │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Species Identification                      │
│  (Mash distance → species_mash_sketches.msh)                    │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    MLST Sequence Typing                          │
│  (BLAST → 7 loci → ST profile matching)                         │
│  Only runs if species is in CdSC                                │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Tox Allele Detection                        │
│  (BLAST → tox allele database → NTTB prediction)                │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│               AMR/Virulence Gene Screening                       │
│  (AMRFinderPlus → custom C. diphtheriae database)               │
│  Includes genomic context analysis                              │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│              Integron Detection (Optional)                       │
│  (Integron_Finder → CALIN, complete, In0 counts)                │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│              Tree Generation (Optional)                          │
│  (JolyTree → Newick tree file)                                  │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Results Aggregation                           │
│  (Combine per-sample results into DataFrame)                    │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Summary Statistics Calculation                  │
│  (summary.py → species counts, tox prevalence, AMR stats)       │
└─────────────────────────────────────────────────────────────────┘
         │
         ├──────────────────┬──────────────────┬──────────────────┐
         ▼                  ▼                  ▼                  ▼
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│   TSV Output    │ │   JSON Output   │ │   HTML Report   │ │   iTOL Files    │
│  {outdir}.txt   │ │  {outdir}.json  │ │  {outdir}.html  │ │ (visualization) │
│  (--format tsv) │ │ (--format json) │ │ (--format html) │ │                 │
└─────────────────┘ └─────────────────┘ └─────────────────┘ └─────────────────┘

Output Files:
  - {outdir}.txt         (TSV results table, default)
  - {outdir}.json        (JSON results, with --format json/all)
  - {outdir}.html        (HTML report, with --format html/all)
  - {strain}.prot.fa     (extracted protein sequences)
  - {strain}.blast.out   (AMRFinderPlus raw output)
  - iTOL template files  (visualization data)
```

## External Tool Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| **Mash** | >= 2.1 | Fast pairwise genome distance estimation |
| **BLAST+** | >= 2.12 | Sequence alignment (blastn, blastp, makeblastdb) |
| **AMRFinderPlus** | >= 3.10 | AMR and virulence gene detection |
| **Integron_Finder** | >= 2.0.5 | Integron detection (optional) |
| **JolyTree** | - | Phylogenetic tree construction (optional) |
| **HMMer** | - | hmmsearch for protein domain searches |
| **Prodigal** | - | Gene prediction (for Integron_Finder) |
| **Infernal** | - | cmsearch for RNA structure searches |

### Tool Availability Check

On startup, diphtOscan validates required tools are available in PATH:
- Core: `mash`, `amrfinder`, `hmmsearch`, `makeblastdb`, `blastn`, `blastp`
- JolyTree: `JolyTree.sh`, `gawk`, `fastme`, `REQ`
- Integron_Finder: `hmmsearch`, `cmsearch`, `prodigal`

## Database Structure

See [DATABASE.md](DATABASE.md) for detailed database documentation.

```
data/
├── mlst/           # MLST alleles and ST profiles
├── tox/            # Tox allele sequences and profiles
├── species/        # Mash sketch reference database
└── resistance/     # AMRFinderPlus + custom C. diphtheriae data
```

## Configuration

diphtOscan uses command-line arguments for all configuration. Key parameters:
- `--min_identity`: Minimum BLAST alignment identity (default: 80%)
- `--min_coverage`: Minimum BLAST alignment coverage (default: 50%)
- `--threads`: Number of parallel processing threads (default: 4)

## Error Handling

- Missing dependencies trigger warnings and graceful degradation
- Invalid input files are reported with descriptive error messages
- Database update failures are caught and reported
- Output directory conflicts can be resolved with `--overwrite`
