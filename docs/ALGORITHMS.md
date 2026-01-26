# diphtOscan Algorithms

This document describes the scientific methodology and algorithms used in diphtOscan.

## Species Identification

### Method: Mash Distance Estimation

diphtOscan uses [Mash](https://github.com/marbl/Mash) for rapid species identification based on MinHash distance estimation.

**Algorithm:**
1. Compute Mash distance between query assembly and reference species sketches
2. Select species with minimum distance
3. Classify match confidence based on distance thresholds

**Thresholds:**
| Distance | Classification | Interpretation |
|----------|----------------|----------------|
| <= 0.05 | Strong match | High confidence species assignment |
| <= 0.10 | Weak match | Possible species, recommend verification |
| > 0.10 | Unknown | Species could not be determined |

**Reference Species:**
- *C. diphtheriae*
- *C. ulcerans*
- *C. pseudotuberculosis*
- *C. belfantii*
- *C. rouxii*
- *C. ramonii*
- *C. silvaticum*

---

## MLST Sequence Typing

### Method: BLAST-based Allele Matching

Multi-Locus Sequence Typing uses 7 housekeeping genes from the PubMLST *C. diphtheriae* scheme.

**Loci:**
`atpA`, `dnaE`, `dnaK`, `fusA`, `leuA`, `odhA`, `rpoB`

**Algorithm:**

1. **Allele Detection**
   - Run BLASTN against allele database for each locus
   - Filter hits by minimum identity and coverage thresholds
   - Select best hit per locus based on BLAST score

2. **Exact vs Imprecise Matches**
   - Exact: 100% identity AND full alignment length
   - Imprecise: Marked with `*` suffix

3. **ST Assignment**
   - Concatenate allele numbers
   - Look up ST in profile database
   - If no exact profile match, find closest ST (locus variant)

4. **Required Exact Matches**
   ```python
   required_exact_matches = len(header) // 2  # At least half of 7 loci = 3
   ```
   An ST is only called if at least 3 loci have exact allele matches.

5. **Locus Variants**
   - `-1LV`: One locus differs from closest ST
   - `-2LV`: Two loci differ
   - etc.

**ST = 0 Conditions:**
- Fewer than 3 exact allele matches
- No close ST variant found
- Allele profile not in database

---

## Tox Gene Detection

### Method: Allele Matching with NTTB Prediction

The tox (diphtheria toxin) gene is detected using BLAST against the PubMLST tox allele database.

**Algorithm:**

1. **Allele Detection**
   - BLAST query against tox allele sequences
   - Apply coverage and identity thresholds
   - Assign allele number

2. **NTTB (Non-Toxigenic Tox-Bearing) Prediction**

   A strain is classified as NTTB when:
   ```
   tox gene detected
   AND coverage < 100%
   AND NOT truncated at contig boundary
   ```

   This indicates a potentially non-functional tox gene due to:
   - Deletions within the gene
   - Frameshift mutations
   - Other disrupting mutations

3. **Truncation at Contig Boundary**

   Genes truncated at contig edges are NOT classified as NTTB because:
   - The full gene may exist but is split across contigs
   - Assembly artifacts should not be interpreted as biological disruption

---

## Truncation Detection

### Method: Amino Acid Coverage Analysis

Truncation is evaluated at the protein level to detect premature stop codons.

**Algorithm:**

1. **Translation**
   - Extract nucleotide sequence from BLAST hit
   - Translate using bacterial genetic code (table 11)
   - Stop at first stop codon

2. **Coverage Calculation**
   ```python
   ref_aa_length = (ref_nucleotide_length - 3) // 3  # Exclude stop codon
   coverage = 100.0 * len(translation) / ref_aa_length
   ```

3. **Truncation Threshold**
   ```python
   cov_threshold = 90.0  # Amino acid coverage threshold (not nucleotide)
   ```
   - Coverage >= 90% at the amino acid level: Gene considered intact
   - Coverage < 90% at the amino acid level: Marked as truncated with percentage

   Note: This threshold applies to the translated protein sequence, not the nucleotide sequence.

4. **Requirements**
   - Hit must start at position 1 of reference
   - Ambiguous bases (non-ACGT) cause sequence truncation at that point

---

## Genomic Context Analysis

### Method: Co-location Detection

Resistance genes are grouped by their proximity on contigs.

**Algorithm:**

1. **Gene Grouping**
   - Group AMR genes by contig
   - Sort genes by start position on each contig

2. **Distance Calculation**
   ```python
   distance = start_position[gene_i+1] - stop_position[gene_i]
   ```

3. **Co-location Threshold**
   ```python
   max_distance = 8000  # base pairs
   ```
   - Genes within 8000 bp are considered co-located (joined by `;`)
   - Genes beyond 8000 bp are separated (joined by `||`)

**Output Format:**
```
geneA;geneB;geneC || geneD
```
Where geneA, geneB, geneC are within 8kb of each other.

---

## AMR/Virulence Gene Detection

### Method: AMRFinderPlus with Custom Database

diphtOscan uses [AMRFinderPlus](https://github.com/ncbi/amr) with an extended database.

**Parameters:**
| Parameter | Value | Description |
|-----------|-------|-------------|
| `--ident_min` | -1 (default) | Uses AMRFinderPlus default identity |
| `--coverage_min` | User-defined (default 50%) | Minimum coverage threshold |
| `--organism` | Corynebacterium_diphtheriae | Species-specific genes |
| `--plus` | Enabled | Include stress, virulence, biocide genes |
| `--translation_table` | 11 | Bacterial genetic code |

**Custom Database Extensions:**
- *C. diphtheriae*-specific virulence factors
- Pili gene clusters (SpaA, SpaD, SpaH types)
- Iron uptake systems
- Additional toxins

---

## Contig Edge Detection

### Method: Position-based Boundary Check

Determines if a partial gene match is due to contig boundary truncation.

**Algorithm:**

```python
def is_contig_edge(hit):
    ref_length = reference_sequence_length * 3  # Convert AA to nucleotides
    seq_length = stop_position - start_position + 1

    if seq_length < ref_length:
        missing_nucleotides = ref_length - seq_length

        # Check if gene extends past contig start
        over_start = (start_position - missing_nucleotides) < 0

        # Check if gene extends past contig end
        over_stop = (contig_length - (stop_position + missing_nucleotides)) < 0

        return over_start or over_stop

    return False
```

Genes truncated at contig edges are annotated with `_end_of_contig`.

---

## BLAST Parameters

### Nucleotide BLAST (blastn)

```bash
blastn -task blastn -db {database} -query {query} \
    -outfmt '6 sacc pident slen length bitscore qseq sstrand sstart send qacc qstart qend qframe' \
    -dust no \
    -evalue 1E-20 \
    -word_size 32 \
    -max_target_seqs 10000 \
    -perc_identity {min_identity}
```

**Key Parameters:**
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `-evalue` | 1E-20 | Stringent threshold for significant hits |
| `-word_size` | 32 | Optimized for similar sequences |
| `-dust no` | Disabled | Don't mask low-complexity regions |

### Hit Filtering

1. **Identity filter**: Remove hits below minimum percent identity
2. **Coverage filter**: Remove hits below minimum reference coverage
3. **Redundancy culling**: Keep only best hit per genomic region
4. **Quality scoring**: `quality = pcid * score * ref_cov`

---

## Integron Detection

### Method: Integron_Finder

When enabled (`-integron`), diphtOscan uses [Integron_Finder](https://github.com/gem-pasteur/Integron_Finder).

**Output Metrics:**
| Type | Description |
|------|-------------|
| `CALIN` | Cluster of Attc sites Lacking Integrase Nearby |
| `complete` | Complete integrons with integrase and attC sites |
| `In0` | Integrase only (no attC cassettes) |

---

## Phylogenetic Tree Generation

### Method: JolyTree

When enabled (`-tree`) with >= 4 assemblies, diphtOscan generates a phylogenetic tree using [JolyTree](https://gitlab.pasteur.fr/GIPhy/JolyTree).

**Requirements:**
- Minimum 4 input assemblies
- Dependencies: JolyTree.sh, gawk, fastme, REQ

**Output:**
- Newick format tree file

---

## References

### Core Tools

1. **Mash** - Ondov BD, Treangen TJ, Melsted P, et al. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*. 2016;17:132. https://doi.org/10.1186/s13059-016-0997-x

2. **AMRFinderPlus** - Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*. 2021;11:12728. https://doi.org/10.1038/s41598-021-91456-0

3. **JolyTree** - Criscuolo A. A fast alignment-free bioinformatics procedure to infer accurate distance-based phylogenetic trees from genome assemblies. *Research Ideas and Outcomes*. 2019;5:e36178. https://doi.org/10.3897/rio.5.e36178

4. **Integron_Finder** - Néron B, Littner E, Haudiquet M, et al. IntegronFinder 2.0: Identification and Analysis of Integrons across Bacteria, with a Focus on Antibiotic Resistance in Klebsiella. *Microorganisms*. 2022;10:700. https://doi.org/10.3390/microorganisms10040700

### diphtOscan Publication

- Hennart M, Panunzi LG, Dupuy B, et al. Population genomics and antimicrobial resistance in *Corynebacterium diphtheriae*. *Peer Community Journal*. 2023;3:e63. https://doi.org/10.24072/pcjournal.307
