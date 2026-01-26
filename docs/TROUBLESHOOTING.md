# Troubleshooting Guide

This document addresses common issues and frequently asked questions about diphtOscan.

## Common Errors

### Dependency Errors

#### "mash missing in path!"
```
/!\ Warning /!\ : mash missing in path!
```

**Solution:**
1. Ensure Mash is installed: `conda install -c bioconda mash`
2. Verify it's in your PATH: `which mash`
3. Make sure the conda environment is activated: `conda activate diphtoscan`

#### "amrfinder missing in path!"
```
/!\ Warning /!\ : amrfinder missing in path!
```

**Solution:**
1. Install AMRFinderPlus: `conda install -c bioconda ncbi-amrfinderplus`
2. Initialize the database: `amrfinder -u` (if using standalone AMRFinderPlus)
3. diphtOscan manages its own database; run `diphtoscan -u` first

#### "Integron_finder missing in path!"
```
/!\ Warning /!\ : Integron_finder missing in path!
```

**Solution:**
1. Integron_Finder is optional; integron analysis will be skipped
2. To enable: `conda install -c bioconda integron_finder>=2.0.5`
3. Ensure version >= 2.0.5 for compatibility

#### "JolyTree.sh missing!"
```
/!\ Warning /!\ : JolyTree.sh missing in /diphtoscan/script/JolyTree/ !
```

**Solution:**
1. Tree generation is optional; analysis continues without it
2. Install JolyTree from: https://gitlab.pasteur.fr/GIPhy/JolyTree
3. Ensure all JolyTree dependencies are installed: gawk, fastme, REQ

### Database Errors

#### "Database not found" or empty results
**Causes:**
- Database not initialized
- Database files corrupted
- Wrong AMRFinderPlus version

**Solution:**
```bash
# Update/initialize the database
diphtoscan -u
```

#### AMRFinderPlus version mismatch
```
Unsupported AMRFinderPlus version: X
```

**Solution:**
- diphtOscan supports AMRFinderPlus versions 3.x and 4.x
- Check version: `amrfinder --version`
- Update if needed: `conda update ncbi-amrfinderplus`

### File/Directory Errors

#### "Directory cannot be created"
```
Directory 'results_xxx' can not be created
```

**Solution:**
1. Check write permissions in the current directory
2. Use `-o` to specify a different output location
3. If directory exists, use `--overwrite` to replace it

#### "Assembly file does not exist"
```
FileNotFoundError: Assembly file X does not exist.
```

**Solution:**
1. Verify the file path is correct
2. Use absolute paths for reliability
3. Check file permissions

### Memory/Performance Issues

#### Analysis is very slow
**Possible causes:**
- Large number of assemblies
- Low thread count
- Large assembly files

**Solution:**
1. Increase threads: `--threads 8`
2. Process assemblies in batches
3. Ensure sufficient RAM (recommended: 8GB+)

#### Out of memory
**Solution:**
1. Reduce number of concurrent assemblies
2. Increase system swap space
3. Process assemblies in smaller batches

---

## Frequently Asked Questions

### Interpreting Results

#### Q: What does `-NTTB` mean?
**A:** NTTB stands for "Non-Toxigenic Tox-Bearing". This annotation indicates:
- The tox gene was detected
- The gene appears to be disrupted (< 100% coverage)
- The gene is NOT truncated at a contig boundary

NTTB strains carry the tox gene but likely cannot produce functional toxin due to mutations, deletions, or frameshifts.

#### Q: What does `*` after a gene name mean?
**A:** The `*` symbol indicates an imprecise match:
- For alleles: Not a 100% identical match to any known allele
- For genes: Found by BLAST rather than exact matching

Example: `ermX*` means ermX was detected but with some sequence differences from the reference.

#### Q: What does `-2LV` mean in ST results?
**A:** The `-nLV` suffix indicates a "Locus Variant":
- `-1LV`: 1 allele differs from the closest known ST
- `-2LV`: 2 alleles differ from the closest known ST
- etc.

Example: `ST8-2LV` means the strain is closest to ST8 but differs at 2 of the 7 MLST loci.

#### Q: Why is ST shown as "0"?
**A:** ST = 0 indicates no sequence type could be assigned. Possible reasons:
1. Fewer than 3 alleles had exact matches (threshold: half of 7 loci)
2. Allele combination not in the database
3. Novel sequence type (consider submitting to PubMLST)

#### Q: Why is ST shown as "NA"?
**A:** "NA" (Not Applicable) means:
- The species is not in the *C. diphtheriae* species complex
- MLST typing only applies to CdSC members

#### Q: What does `_end_of_contig` mean?
**A:** This annotation indicates the gene appears truncated because it reaches the end of a contig. This could mean:
- The gene is split across contigs (assembly fragmentation)
- The gene is actually truncated at this position
- Further investigation or better assembly may be needed

#### Q: How should I interpret genomic context?
**A:** The genomic context shows spatial relationships between resistance genes:
```
ermX;aph(3')-Ia || tetA
```
- `ermX` and `aph(3')-Ia` are close together (within 8000 bp)
- `tetA` is on a different contig or far away (> 8000 bp)

Co-located genes may be on the same mobile element.

### Analysis Parameters

#### Q: What are recommended settings for high-quality assemblies?
```bash
diphtoscan -a *.fasta --mlst --tox --resistance_virulence \
    --min_identity 90 --min_coverage 80 --threads 8
```

#### Q: What settings should I use for draft assemblies?
```bash
diphtoscan -a *.fasta --mlst --tox --resistance_virulence \
    --min_identity 80 --min_coverage 50 --threads 8
```
Lower thresholds accommodate fragmented assemblies but may increase false positives.

#### Q: How do I analyze many genomes efficiently?
```bash
# Run with maximum threads
diphtoscan -a genomes/*.fasta --mlst --resistance_virulence \
    --threads $(nproc) -o batch_results
```

### Database Management

#### Q: How often should I update the database?
**A:** Update before important analyses or at least monthly:
```bash
diphtoscan -u
```
This fetches the latest MLST profiles and AMRFinderPlus data.

#### Q: Can I use a custom database?
**A:** The custom *C. diphtheriae* database is automatically merged with AMRFinderPlus. For additional customizations, modify files in `data/resistance/Corynebacterium_diphtheriae/`.

---

## Assembly Quality Requirements

### Minimum Requirements

| Metric | Minimum | Recommended |
|--------|---------|-------------|
| N50 | 10,000 bp | > 50,000 bp |
| Total length | 2.0 Mb | 2.3-2.6 Mb |
| Contigs | < 500 | < 100 |
| Coverage | 20x | > 50x |

### Quality Impact on Results

| Issue | Impact |
|-------|--------|
| Low N50 | More `_end_of_contig` annotations, incomplete MLST |
| High contamination | False positive gene detection, wrong species |
| Low coverage | Missing genes, partial matches |
| Adapter contamination | BLAST artifacts, false positives |

### Quality Check

Before running diphtOscan, consider:
```bash
# Check assembly stats with seqkit
seqkit stats assembly.fasta

# Check for contamination with CheckM or similar tools
```

---

## Docker-Specific Issues

### Container exits immediately

**Cause:** Running without arguments shows help and exits.

**Solution:**
```bash
# Correct: provide analysis arguments
docker run --rm -v /data:/data diphtoscan -a /data/assembly.fasta --mlst
```

### Permission denied on output files

**Cause:** Container user (UID 1000) doesn't match host user.

**Solution:**
```bash
# Option 1: Run as your user
docker run --rm --user $(id -u):$(id -g) -v /data:/data diphtoscan ...

# Option 2: Fix permissions after run
sudo chown -R $(id -u):$(id -g) /path/to/output
```

### Database not persisted between runs

**Cause:** Database is inside the container and lost when container stops.

**Solution:**
```bash
# Use a named volume for the database
docker run --rm -v diphtoscan-db:/opt/diphtoscan/data diphtoscan -u

# Use the same volume when running analyses
docker run --rm -v diphtoscan-db:/opt/diphtoscan/data -v /data:/data diphtoscan -a /data/assembly.fasta --mlst
```

### Out of memory in Docker

**Solution:**
```bash
# Increase Docker memory limit
docker run --rm --memory=8g -v /data:/data diphtoscan -a /data/*.fasta --mlst
```

### Cannot access files in mounted volumes

**Cause:** Paths inside container differ from host paths.

**Solution:**
```bash
# Files in /host/path are accessible at /data inside container
docker run --rm -v /host/path:/data diphtoscan -a /data/assembly.fasta ...
#                                                  ^^^^^^ use container path
```

---

## Getting Help

1. **Check this guide** for common issues
2. **Review the documentation** in the `docs/` folder
3. **GitHub Issues**: https://github.com/AMR-Hackathon-2025/diphtOscan/issues
4. **Contact**: See README.md for author contact information
