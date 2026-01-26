---
title: 'diphtOscan: Comprehensive genomic analysis of the Corynebacterium diphtheriae species complex'
tags:
  - Python
  - bioinformatics
  - genomics
  - Corynebacterium diphtheriae
  - MLST
  - antimicrobial resistance
authors:
  - name: Melanie Hennart
    affiliation: 1
  - name: Martin Rethoret-Pasty
    affiliation: 1
affiliations:
  - name: Institut Pasteur, Biodiversity and Epidemiology of Bacterial Pathogens, Paris, France
    index: 1
date: January 2026
bibliography: paper.bib
---

# Summary

diphtOscan is a Python command-line tool for comprehensive genomic characterization of *Corynebacterium diphtheriae* and related species within the *C. diphtheriae* species complex (CdSC). The tool analyzes genome assemblies to determine species identity, multi-locus sequence type (MLST), diphtheria toxin gene (*tox*) presence and functionality, antimicrobial resistance determinants, and virulence factors. diphtOscan provides researchers and public health laboratories with a standardized, reproducible workflow for characterizing clinical and environmental CdSC isolates. The tool produces multiple output formats including tab-separated tables, structured JSON for programmatic access, self-contained interactive HTML reports for data exploration, and visualization-ready files compatible with the Interactive Tree of Life (iTOL) platform.

# Statement of Need

Diphtheria, caused primarily by toxigenic strains of *Corynebacterium diphtheriae*, remains a significant global health concern despite widespread vaccine availability [@Hennart2023]. The disease causes respiratory illness characterized by pseudomembrane formation and can lead to systemic toxicity affecting the heart and nervous system. Related species within the CdSC, including *C. ulcerans*, *C. pseudotuberculosis*, *C. belfantii*, *C. rouxii*, *C. ramonii*, and *C. silvaticum*, can also cause diphtheria-like illness in humans and animals, with *C. ulcerans* increasingly recognized as a source of zoonotic infections.

The emergence of antimicrobial-resistant strains poses challenges for treatment, particularly in outbreak settings where empirical therapy may be compromised. Additionally, the identification of non-toxigenic *tox*-bearing (NTTB) isolates, which carry disrupted toxin genes that could potentially revert to functional forms, complicates surveillance efforts and risk assessment.

Whole-genome sequencing has become an essential tool for epidemiological surveillance and outbreak investigation. However, researchers analyzing CdSC genomes face a fragmented landscape of tools: general-purpose software for species identification, separate databases for MLST typing, and resistance gene detection tools that lack CdSC-specific annotations. This fragmentation creates barriers to reproducible analysis and requires substantial bioinformatics expertise to integrate results from multiple sources.

diphtOscan addresses this gap by providing an integrated, species-complex-specific pipeline that automates comprehensive genomic characterization. The tool incorporates curated databases for CdSC virulence factors and produces standardized output suitable for both individual isolate analysis and large-scale surveillance studies. Public health laboratories can use diphtOscan to rapidly characterize clinical isolates, while researchers can process hundreds of genomes for population-level studies.

# Key Features

## Species Identification

diphtOscan uses Mash [@Ondov2016] for rapid species identification based on MinHash distance estimation against reference genome sketches. Species assignments are classified by confidence level: strong matches (distance <= 0.05), weak matches (distance <= 0.10), or unknown. The tool recognizes all seven currently described species within the CdSC.

## MLST Sequence Typing

Multi-locus sequence typing follows the PubMLST *C. diphtheriae* scheme [@Jolley2010], which employs seven housekeeping genes (*atpA*, *dnaE*, *dnaK*, *fusA*, *leuA*, *odhA*, *rpoB*). diphtOscan performs BLAST-based allele matching, identifies exact and imprecise allele matches, and reports locus variants when the allelic profile differs from known sequence types.

## Diphtheria Toxin Analysis

The *tox* gene is detected using BLAST searches against the PubMLST tox allele database. diphtOscan predicts NTTB status when a *tox* gene is present but incomplete, distinguishing biological disruption (deletions, frameshifts) from assembly artifacts at contig boundaries. Protein-level truncation analysis translates sequences and calculates amino acid coverage to identify premature stop codons.

## Antimicrobial Resistance and Virulence Screening

diphtOscan integrates AMRFinderPlus [@Feldgarden2021] with extended annotations for CdSC-specific virulence factors. The custom database includes pili gene clusters (SpaA, SpaD, SpaH types), iron uptake systems, biovar-associated genes (*spuA*, nitrate reductase cluster), and additional toxins beyond the core AMRFinderPlus catalog.

## Genomic Context Analysis

Resistance and virulence genes are analyzed for genomic co-location. Genes within 8,000 base pairs on the same contig are reported as potentially co-transferred elements, facilitating identification of resistance islands and mobile genetic elements.

## Assembly Quality Control

diphtOscan calculates assembly quality metrics for input genomes, including total assembly length, contig count, N50 and L50 statistics, GC content, and ambiguous base percentage. These metrics help users identify poor-quality assemblies that may produce unreliable results. Quality control can be disabled with the `--no-qc` flag for pre-validated assemblies.

## Output Formats and Reporting

diphtOscan supports multiple output formats to accommodate different user needs:

- **TSV format**: Tab-separated tables suitable for spreadsheet analysis and downstream processing
- **JSON format**: Structured data with metadata and summary statistics for programmatic access and integration with analysis pipelines
- **HTML reports**: Self-contained interactive reports featuring summary cards, searchable and sortable results tables, color-coded gene presence visualization, and species distribution charts

The HTML reports require no external dependencies and can be shared as single files, making them ideal for collaborator communication and archival purposes. Users can generate all formats simultaneously with the `--format all` option.

## Batch Processing and User Experience

For large-scale analyses, diphtOscan displays progress indicators when processing multiple assemblies and generates summary statistics showing species distribution, toxin gene prevalence, and antimicrobial resistance patterns across the sample set. Verbose (`-v`) and quiet (`-q`) modes allow users to control output verbosity according to their workflow requirements.

## Visualization Output

diphtOscan generates iTOL-compatible annotation files for phylogenetic tree visualization, including binary presence/absence datasets, color strip annotations for sequence types, and AMR gene family summaries. When four or more assemblies are provided, optional phylogenetic tree construction is available via JolyTree [@Criscuolo2019].

## Database Updates

The tool provides automated database updates from PubMLST for MLST alleles and profiles, and from NCBI for AMRFinderPlus data, ensuring analyses use current reference information.

# Implementation

diphtOscan is implemented in Python 3.8+ and orchestrates established bioinformatics tools including Mash for species identification, BLAST+ for sequence alignment, and AMRFinderPlus for resistance gene detection. Optional modules include Integron_Finder [@Neron2022] for integron detection and JolyTree for phylogenetic tree construction. The modular architecture allows users to select specific analyses based on their needs, from quick species identification to full characterization with resistance profiling.

The command-line interface employs a subcommand structure (`diphtoscan species`, `diphtoscan mlst`, `diphtoscan amr`, `diphtoscan tox`, `diphtoscan qc`, `diphtoscan all`) that enables focused analyses while maintaining backward compatibility with earlier versions. This design reduces execution time and dependency requirements when only specific analyses are needed.

The tool accepts one or more genome assemblies in FASTA format and produces results in multiple formats. Tab-separated tables provide compatibility with spreadsheet software and Unix text processing tools. JSON output includes structured metadata and summary statistics suitable for integration with laboratory information management systems and automated reporting pipelines. Interactive HTML reports, generated using Jinja2 templates, provide immediate data exploration capabilities without requiring additional software.

Installation is supported via conda for dependency management, with the core tool installed via pip. A Docker image is available for containerized execution in computing clusters and cloud environments. Configuration is handled entirely through command-line arguments, supporting reproducible analysis pipelines and integration into automated workflows. Users can adjust minimum identity and coverage thresholds for BLAST searches and specify the number of threads for parallel processing. A comprehensive test suite with over 330 unit and integration tests ensures software reliability across updates.

# Mentions of Related Work

diphtOscan's architecture draws inspiration from Kleborate, a similar species-specific typing tool for *Klebsiella pneumoniae*. While general-purpose tools like mlst and AMRFinderPlus provide individual analysis components, diphtOscan integrates these into a CdSC-specific workflow with custom databases and NTTB prediction capabilities not available elsewhere. The tool complements PubMLST web-based typing by providing programmatic access suitable for batch processing and integration into bioinformatics pipelines.

# Acknowledgements

We thank the Biodiversity and Epidemiology of Bacterial Pathogens team at Institut Pasteur for testing and feedback. We acknowledge the PubMLST database and NCBI AMRFinderPlus teams for maintaining the reference databases that diphtOscan relies upon.

# References
