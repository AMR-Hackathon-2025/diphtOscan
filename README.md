## _diphtOscan_

_diphtOscan_ is a command line script written in [Python](https://www.python.org/). _DiphtOscan_ runs on UNIX, Linux and most OS X operating systems.
For more details, see the associated publication (https://www.biorxiv.org/content/10.1101/2023.02.20.529124v1xxx).

_diphtOscan_ is a tool to search genomic assemblies of _Corynebacterium diphtheriae_ and other species of the _Corynebacterium diphtheriae_ species complex (CdSC) for:
* Species (e.g. _C. diphtheriae_, _C. ulcerans_, _C. pseudotuberculosis_, _C. belfantii_, _C. rouxii_ , _C. ramonii_ and _C. silvaticum_)
* Biovar-associated genes (_spuA_, nitrate reductase gene cluster)
* MLST sequence type
* Virulence factors, including _tox_ gene detection and disruption prediction
* Antimicrobial resistance determinants: acquired genes (_ermX_, _pbp2m_, …) and SNPs (e.g., _rpoB_, _gyrA_)
* Genomic context of genomic features associated with resistance
* Presence of integrons (using Integron Finder: https://github.com/gem-pasteur/Integron_Finder) 
* Tree building (using JolyTree: https://gitlab.pasteur.fr/GIPhy/JolyTree)

## Installation and execution

**A.** Install the following programs and tools, or verify that they are already installed with the required version:
* [python](https://www.python.org/) version >= 3.8.3
* [mash](http://mash.readthedocs.io/en/latest/) version >= 2.3 
* [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version >= 2.13.0
* [amrfinder](https://github.com/ncbi/amr/wiki) version >= 3.12.8
* [hmmer](http://hmmer.org/)

Note: 
Need for the `--tree` i.e. using the Jolytree tool (https://gitlab.pasteur.fr/GIPhy/JolyTree)
* [Jolytree](https://gitlab.pasteur.fr/GIPhy/JolyTree) version >= 2.1 (The JollyTree repository should be cloned into diphtoscan/script/)

Need for the `--integron` i.e. using the Integron Finder version 2 tool (https://github.com/gem-pasteur/Integron_Finder)

* [integron_finder](https://github.com/gem-pasteur/Integron_Finder) version >= 2.0.5


**B.** Clone this repository with the following command line:
```bash
git clone https://gitlab.pasteur.fr/BEBP/diphtoscan.git
```

**C.** Give the execute permission to the file `diphtoscan/__main__.py`:
```bash
chmod +x diphtoscan/__main__.py
```

**D.** Run the installation of the JollyTree module (if you plan to use it) :
```bash
diphtoscan/script/update_tools.sh
```

**E.** Execute _diphtOscan_ with the following command line model:
```bash
python __main__.py [options]
```

## Usage

Launch _diphtOscan_ without option to read the following documentation:

```
usage: __main__.py -a ASSEMBLIES [ASSEMBLIES ...] [-u] [-st] [-t] [-res_vir] [-plus] [-integron] [-o OUTDIR]
                   [--min_identity MIN_IDENTITY] [--min_coverage MIN_COVERAGE] [--threads THREADS] [-tree] 
                   [--overwrite] [-h] [--version]

diphtOscan: a tool for characterising virulence and resistance in Corynebacterium

Updating option:
  -u, --update          Update database MLST, Tox Allele & AMR (default: no).The database update can be executed on its own without the -a option.

Required arguments:
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        FASTA file(s) for assemblies.

Screening options:
  -st, --mlst           Turn on species Corynebacterium diphtheriae species complex (CdSC) and MLST sequence type
                        (default: no)
  -t, --tox             Turn on tox allele (default: no)
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
  --overwrite           Allows the output directory to be overwritten if it already exists

Phylogenetic tree:
  -tree, --tree         Generates a phylogenetic tree from JolyTree

Help:
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
```

## Example

In order to illustrate the usefulness of _diphtOscan_ and to describe its output files, the following use case example describes its usage for inferring a phylogenetic tree of _Corynebacterium diphtheriae_ genomes derived from the analysis of [Hennart et al](https://peercommunityjournal.org/articles/10.24072/pcjournal.307/).

##### Running _diphtOscan_

The following command line allows the script `diphtOscan` to be launched with default options on 8 threads:
```bash
python3 __main__.py -a $genomes --mlst --resistance_virulence --threads 8 -o Cdiphteriae
```

As the basename was set to 'Cdiphteriae', _diphtOscan_ writes in few minutes the four following output files:

* `Cdiphteriae.csv`: result file 
* `$strain.fa`: extracted sequences (for every assemblie files) 
* `$strain.out`: BLAST output file (for every assemblie files) 





