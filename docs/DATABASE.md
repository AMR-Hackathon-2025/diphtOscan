# diphtOscan Database Documentation

This document describes the database structure, update procedures, and data sources used by diphtOscan.

## Directory Structure

```
data/
├── mlst/                    # MLST scheme data
│   ├── sequences/           # Individual allele FASTA files
│   │   ├── atpA.fas
│   │   ├── dnaE.fas
│   │   ├── dnaK.fas
│   │   ├── fusA.fas
│   │   ├── leuA.fas
│   │   ├── odhA.fas
│   │   └── rpoB.fas
│   ├── pubmlst_diphtheria_seqdef_scheme_3.fas  # Combined allele database
│   └── st_profiles.txt      # ST profile definitions
│
├── tox/                     # Tox allele scheme data
│   ├── sequences/           # Tox allele FASTA files
│   │   └── tox.fas
│   ├── pubmlst_diphtheria_seqdef_scheme_4.fas  # Tox allele database
│   └── tox_profiles.txt     # Tox profile definitions
│
├── species/                 # Species identification data
│   └── species_mash_sketches.msh  # Mash sketch database
│
└── resistance/              # AMR and virulence data
    ├── Corynebacterium_diphtheriae/  # Custom C. diphtheriae data
    │   ├── AMRProt_Cd       # Custom protein sequences
    │   └── fam_Cd.tab       # Custom gene family annotations
    │
    └── YYYY-MM-DD/          # AMRFinderPlus database (dated)
        ├── AMRProt          # Protein sequences (v3) or AMRProt.fa (v4)
        ├── AMRProt.p*       # BLAST protein database files
        ├── fam.tab/tsv      # Gene family annotations
        ├── version.txt      # Database version
        └── ...              # Other AMRFinderPlus files
```

## Data Sources

### MLST Data (PubMLST)

**Source:** [PubMLST *C. diphtheriae* database](https://pubmlst.org/organisms/corynebacterium-diphtheriae)

**API Endpoint:** `https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef`

**Scheme 3 - MLST:**
- 7 housekeeping gene loci
- Sequence type profiles
- Updated via `diphtoscan -u`

**Scheme 4 - Tox:**
- Tox allele sequences
- Tox profiles
- Updated via `diphtoscan -u`

### Species Reference Data

**Source:** Pre-computed Mash sketches of reference genomes

**Species included:**
- *Corynebacterium diphtheriae*
- *Corynebacterium ulcerans*
- *Corynebacterium pseudotuberculosis*
- *Corynebacterium belfantii*
- *Corynebacterium rouxii*
- *Corynebacterium ramonii*
- *Corynebacterium silvaticum*

**Format:** Mash sketch file (`.msh`)

### AMRFinderPlus Database

**Source:** [NCBI AMRFinderPlus](https://github.com/ncbi/amr)

**URLs:**
- AMRFinderPlus v3: `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/`
- AMRFinderPlus v4: `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/`

**Contents:**
- AMR gene protein sequences
- Virulence factor sequences
- Gene family classifications
- Point mutation data

### Custom *C. diphtheriae* Data

**Source:** Institut Pasteur BEBP team

**Contents:**
Custom annotations for *C. diphtheriae*-specific:
- Virulence factors (pili, adhesins)
- Iron uptake systems
- Toxin variants
- Species-specific resistance genes

## Database Update Process

### Automatic Update

```bash
diphtoscan -u
```

This command performs the following steps:

1. **MLST Database Update**
   - Removes old MLST files
   - Downloads allele sequences from PubMLST API
   - Downloads ST profiles
   - Creates combined FASTA database

2. **Tox Database Update**
   - Removes old tox files
   - Downloads tox allele sequences
   - Downloads tox profiles

3. **AMRFinderPlus Database Update**
   - Downloads latest AMRFinderPlus database
   - Appends custom *C. diphtheriae* data
   - Builds BLAST protein database

### Manual Update Steps

If automatic update fails, you can manually update components:

```bash
# 1. Update MLST (requires internet access to PubMLST)
# The update function handles this automatically

# 2. Update AMRFinderPlus database manually
amrfinder -u  # Updates the system AMRFinderPlus database

# 3. Then run diphtoscan update to merge custom data
diphtoscan -u
```

## AMRFinderPlus Version Compatibility

### Version 3.x

- Database format: `.tab` suffix for TSV files
- Protein file: `AMRProt` (no extension)
- Column names: `Gene symbol`, `% Coverage of reference sequence`

### Version 4.x

- Database format: `.tsv` suffix for TSV files
- Protein file: `AMRProt.fa`
- Column names: `Element symbol`, `% Coverage of reference`

diphtOscan automatically detects the AMRFinderPlus version and uses appropriate file names and column names.

## Custom Gene Classifications

The following gene families are classified in the custom database:

### Virulence Factors

| Category | Genes |
|----------|-------|
| TOXIN | tox |
| OTHER_TOXINS | pld (phospholipase D) |
| SpaA-type pili | spaA, spaB, spaC, srtA |
| SpaD-type pili | spaD, spaE, spaF, srtB, srtC |
| SpaH-type pili | spaG, spaH, spaI, srtD, srtE |
| VIRULENCE/ADHESIN | cbpA, nanH |

### Iron Uptake Systems

- irp1ABCD
- irp2ABCDEFGHI
- irp2JKLMN
- irp6ABC
- iusABCDE
- iutABCDE
- htaA-hmuTUV-htaBC
- hmuO
- frgCBAD
- ciuABCD
- ciuEFG
- chtAB
- chtC
- cdtQP-sidBA-ddpABCD
- HbpA

### Biovar Markers

- spuA (spermidine/putrescine catabolism)
- narG, narIJHK (nitrate reductase cluster)

## Database File Formats

### ST Profiles (`st_profiles.txt`)

Tab-separated file:
```
ST	atpA	dnaE	dnaK	fusA	leuA	odhA	rpoB
1	1	1	1	1	1	1	1
2	1	2	1	1	1	1	1
...
```

### Allele FASTA Files

Standard FASTA format with locus-specific headers:
```
>atpA_1
ATGAAACGTGCAGTTCGTGAA...
>atpA_2
ATGAAACGTGCAGTTCGTGAG...
```

### Family Annotations (`fam.tab`)

Tab-separated file with gene metadata:
```
#node_id	parent_node_id	class	subclass	...
tox	VIRULENCE_Cdiphth	TOXIN	TOXIN	...
spaA	VIRULENCE_Cdiphth	SpaA-type_pili_diphtheriae	...
```

## Troubleshooting Database Issues

### Database Not Found

```bash
# Verify database directory exists
ls -la /path/to/diphtoscan/data/

# Re-run update
diphtoscan -u
```

### Corrupted Database

```bash
# Remove and regenerate
rm -rf data/mlst/sequences data/mlst/*.fas data/mlst/st_profiles.txt
rm -rf data/tox/sequences data/tox/*.fas data/tox/tox_profiles.txt
rm -rf data/resistance/202*  # Remove dated AMR databases
diphtoscan -u
```

### Network Errors During Update

If update fails due to network issues:
1. Check internet connectivity
2. Verify access to `bigsdb.pasteur.fr` and `ftp.ncbi.nlm.nih.gov`
3. Retry the update

### Version Mismatch

If you see errors about missing columns or files:
1. Check AMRFinderPlus version: `amrfinder --version`
2. Ensure version is 3.x or 4.x
3. Re-run database update: `diphtoscan -u`
