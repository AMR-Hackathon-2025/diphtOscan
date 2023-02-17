# diphtOscan

_diphtOscan_ is a command line script written in [Python](https://www.python.org/). _DiphtOscan_ runs on UNIX, Linux and most OS X operating systems.
For more details, see the associated publication (xxx).

_diphtOscan_  is a tool to screen genome assemblies of _Corynebacterium diphtheriae_ and the _Corynebacterium diphtheriae_ species complex (Cdc) for:
 * Species (e.g. _C. diphtheriae_, _C. belfantii_, _C. rouxii_, _C. ulcerans_, _C. silvaticum_ and _C. pseudotuberculosis_)
 * MLST sequence type
 * Virulence loci 
 * Antimicrobial resistance: acquired genes, SNPs
 * Biovar prediction
 * Detection of tox gene (presence/absence and get _tox_ alleles)
 * Get the genomic context of the Antimicrobial resistance genes


## Installation and execution

**A.** Install the following programs and tools, or verify that they are already installed with the required version:
* [python](https://www.python.org/) version >= 3.8.3
* [mash](http://mash.readthedocs.io/en/latest/) version >= 2.3; 
* [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version >= 2.13.0
* [hmmer](http://hmmer.org/download.html) version >= 3.3.2
* [amrfinder](https://github.com/ncbi/amr/wiki) version >= 3.11.2

Note: 
Need for the `--tree` i.e. using the Jolytree tool (https://gitlab.pasteur.fr/GIPhy/JolyTree)
* [fastme] version >= 2.1.6.4 
* [REQ] version >= 1.3
* [Jolytree](https://gitlab.pasteur.fr/GIPhy/JolyTree) version >= 2.1

Need for the `--integron` i.e. using the Integron Finder version 2 tool (https://github.com/gem-pasteur/Integron_Finder)

* [infernal] version >=
* [prodigal] version >=
* [integron_finder](https://github.com/gem-pasteur/Integron_Finder) version >= 2.0.2


**B.** Clone this repository with the following command line:
```bash
git clone https://gitlab.pasteur.fr/BEBP/diphtoscan.git
```

**C.** Give the execute permission to the file `diphtOscan/__main__.py`:
```bash
chmod +x __main__.py
```

**D.** Execute _JolyTree_ with the following command line model:
```bash
python __main__.py [options]
```

## Usage

Launch _diphtOscan_ without option to read the following documentation:

```
usage: __main__.py -a ASSEMBLIES [ASSEMBLIES ...] [-u] [-t] [-res_vir] [-plus] [-integron] [-o OUTDIR]
                   [--min_identity MIN_IDENTITY] [--min_coverage MIN_COVERAGE] [--threads THREADS] [-tree] [-h]
                   [--version]

diphtOscan: a tool for characterising virulence and resistance in Corynebacterium

Required arguments:
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies

Screening options:
  -u, --update          Update database MLST et AMR (default: no)
  -t, --taxonomy        Turn on species Corynebacterium diphtheriae species complex (CdSC) and MLST sequence type
                        (default: no)
  -res_vir, --resistance_virulence
                        Turn on resistance and virulence genes screening (default: no resistance and virulence gene
                        screening)
  -plus, --extend_genotyping
                        Turn on all virulence genes screening (default: no all virulence gene screening)
  -integron, --integron
                        Screening the intregon(default: no)

Output options:
  -o OUTDIR, --outdir OUTDIR
                        Folder for detailed output (default: results_YYYY-MM-DD_II-MM-SS_PP)

Settings:
  --min_identity MIN_IDENTITY
                        Minimum alignment identity for main results (default: 80)
  --min_coverage MIN_COVERAGE
                        Minimum alignment coverage for main results (default: 50)
  --threads THREADS     The number of threads to use for processing. (default: 4)

Phylogenetic tree:
  -tree, --tree         Generates a phylogenetic tree from JolyTree

Help:
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
```

## Notes


## Example

In order to illustrate the usefulness of _diphtOscan_ and to describe its output files, the following use case example describes its usage for inferring a phylogenetic tree of _Klebsiella_ genomes derived from the analysis of [Rodrigues et al. (2019)](https://doi.org/10.1016/j.resmic.2019.02.003).

##### Running _diphtOscan_

The following command line allows the script `diphtOscan` to be launched with default options on 8 threads:
```bash
python3 __main__.py -a $genomes --taxonomy --resistance_virulence --threads 8 -o Cdiphteriae
```

As the basename was set to 'Cdiphteriae', _diphtOscan_ writes in few minutes the four following output files:

* `Cdiphteriae.csv`: result file 
* `$strain.fa`: extracted sequences (for every assemblie files) 
* `$strain.out`: BLAST output file (for every assemblie files) 



