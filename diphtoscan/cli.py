#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# Copyright (C) 2020  Melanie HENNART                                         #
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.      #
#                                                                             #
#                                                                             #
#  Contact:                                                                   #
#                                                                             #
#    Melanie HENNART, PhD Student                                             #
#    melanie.hennart@pasteur.fr                                               #
#    Biodiversity and Epidemiology of Bacterial Pathogens                     #
#    Institut Pasteur                                                         #
#    25-28, Rue du Docteur Roux                                               #
#    75015 Paris Cedex 15                                                     #
#    France                                                                   #
#                                                                             #
###############################################################################

"""
diphtOscan is a tool to screen genome assemblies of the diphtheriae species
complex (DiphSC) for:
     - Species (e.g. C. diphtheriae, C. belfantii, C. rouxii, C. ulcerans,
                C. ramonii and C. pseudotuberculosis)
     - Biovar-associated genes
     - MLST sequence type
     - Virulence factors
     - Antimicrobial resistance: acquired genes, SNPs & genomic context
     - Detection of tox gene (Presence/Absence & Allele)
     - Detection of integrons (Integron_Finder)
     - Tree building (JolyTree)

Usage:
======
    diphtoscan all -a genome1.fasta genome2.fasta
    diphtoscan species -a genome.fasta
    diphtoscan mlst -a genome.fasta
    diphtoscan amr -a genome.fasta
    diphtoscan tox -a genome.fasta
    diphtoscan qc -a genome.fasta
    diphtoscan update
"""

__authors__ = ("Melanie HENNART; Martin RETHORET-PASTY")
__contact__ = ("martin.rethoret-pasty@pasteur.fr")
__version__ = "1.8.0"
__copyright__ = "copyleft"
__date__ = "2024/03/04"

###############################################################################
#                                                                             #
# ================                                                            #
# = INSTALLATION =                                                            #
# ================                                                            #
#                                                                             #
# [1] REQUIREMENTS =========================================================  #
#                                                                             #
# -- Mash: fast pairwise p-distance estimation -----------------------------  #
#    VERSION >= 2.1                                                           #
#    src: github.com/marbl/Mash                                               #
#    Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S,      #
#      Phillippy AM (2016) Mash: fast  genome  and  metagenome  distance      #
#      estimation using MinHash. Genome Biology, 17:132.                      #
#      doi:10.1186/s13059-016-0997-x                                          #
#                                                                             #
#PATH_MASH="mash"                                                             #
#                                                                             #
#                                                                             #
# [2] NOTES ON THE USE OF JOLYTREE WITH SLURM (slurm.schedmd.com) ==========  #
#                                                                             #
#os.system("module purge") ;
#os.system("module load Mash/2.1") ;                                          #
#PATH_MASH = "module load Mash/2.1"                                           #
#                                                                             #
#                                                                             #
###############################################################################

###############################################################################
# ================                                                            #
# = NOTE         =                                                            #
# ================                                                            #
#                                                                             #
# mash sketch -o reference genome1.fna genome2.fna                            #
# mash info reference.msh                                                     #
#                                                                             #
###############################################################################

import sys
import os
import glob
import json
import pandas as pd
import datetime
import os.path
import argparse
import subprocess
import shutil

from typing import List
from tqdm import tqdm
from .species import get_species_results, is_cd_complex
from .template_iTOL import spuA, narG, toxin, amr_families
from .updating_database import update_database
from .jolytree_generation import generate_jolytree
from .assembly_qc import calculate_assembly_stats, format_qc_results

from .utils import (
    get_chromosome_mlst_results,
    get_tox_results,
    get_chromosome_mlst_header,
    get_tox_header,
    delete_virulence_extended,
    is_non_zero_file,
    armfinder_to_table,
    get_genomic_context,
    find_resistance_db,
    determine_biovar
)
from .logging_config import setup_logging, get_logger, log_info, log_debug, log_warning, log_error
from .summary import calculate_summary_statistics, format_summary_for_console
from .html_report import generate_html_report


def test_unique_dependency(name: str) -> bool:
    """
    Check if a single command-line dependency is available in the system PATH.

    Parameters
    ----------
    name : str
        Name of the executable to check.

    Returns
    -------
    bool
        True if the dependency is found in PATH, False otherwise.
    """
    return shutil.which(name) is not None


def test_multiple_dependencies(dependencies: List[str]) -> None:
    """
    Check if multiple command-line dependencies are available in the system PATH.

    Exits the program with an error if any dependency is missing.

    Parameters
    ----------
    dependencies : List[str]
        List of executable names to check.

    Raises
    ------
    SystemExit
        If any dependency is not found in PATH.
    """
    for dependency in dependencies:
        presence = test_unique_dependency(dependency)
        if presence is not True:
            print(f'/!\\ Warning /!\\ : {dependency} missing in path!')
            sys.exit(-1)


def test_required_dependency(args):
    """
    Validate that all required external dependencies are available.

    Checks for core diphtOscan dependencies (Mash, BLAST+, AMRFinderPlus),
    and optionally checks for Integron_Finder and JolyTree dependencies
    based on command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing flags for optional analyses.

    Returns
    -------
    argparse.Namespace
        The arguments object, potentially with integron and tree flags
        set to False if their dependencies are missing.

    Notes
    -----
    Core dependencies: mash, amrfinder, hmmsearch, makeblastdb, blastn, blastp
    JolyTree dependencies: JolyTree.sh, gawk, fastme, REQ
    Integron_Finder dependencies: hmmsearch, cmsearch, prodigal
    """
    diphtoscan_dependencies = ["mash",'amrfinder','hmmsearch', 'makeblastdb','blastn', 'blastp']
    joly_tree_dependencies = ["JolyTree.sh", "gawk",'fastme','REQ']
    integron_fender_dependencies = ['hmmsearch', 'cmsearch', 'prodigal']

    if args.assemblies is None:  # TODO: Ensure that dependencies are not required to update the database
        return args

    subprocess.run(["echo", "Dependency testing"])
    test_multiple_dependencies(diphtoscan_dependencies)

    if hasattr(args, 'integron') and args.integron:
        rc = test_unique_dependency("Integron_finder")
        test_multiple_dependencies(integron_fender_dependencies)
        if rc == 0:
            args.integron = True
        else:
            print('/!\\ Warning /!\\ : Integron_finder missing in path! Integron analysis not carried out.')
            args.integron = False

    if hasattr(args, 'tree') and args.tree:
        test_multiple_dependencies(joly_tree_dependencies)
        args.tree = True
    else:
        print('/!\\ Warning /!\\ : JolyTree.sh missing in /diphtoscan/script/JolyTree/ ! Joly_tree representation not carried out.')
        if hasattr(args, 'tree'):
            args.tree = False
    print('\n')
    return args


def redefine_output_file(args):
    """
    Set up output directory for overwrite mode.

    When --overwrite is used and the output directory exists, this function
    creates a temporary directory for processing and prepares for final
    file transfer.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments with outdir attribute.

    Returns
    -------
    tuple
        (modified_args, final_output_path) where modified_args has outdir
        set to a temporary directory name.
    """
    if not os.path.exists(args.outdir):
        print(f"Directory {args.outdir} does not exist.")
        try:
            os.makedirs(args.outdir)
            print("Directory '%s' created successfully \n" %args.outdir)
        except OSError :
            print("Directory '%s' can not be created \n"  %args.outdir)
            sys.exit(0)
    final_output_path = args.outdir
    args.outdir = f"{args.outdir}_temp_folder"
    return args, final_output_path


def rename_temp_folder_file(outdir: str, directory: str) -> None:
    """
    Rename the main results file from temporary to final naming convention.

    Removes the '_temp_folder' suffix from the results file name.

    Parameters
    ----------
    outdir : str
        Final output directory path.
    directory : str
        Temporary directory path containing the file to rename.
    """
    path, file_name = os.path.split(directory)
    expected_filename = f"{file_name}.txt"
    file_path = os.path.join(directory, expected_filename)

    print(f"Renaming {file_path} to {expected_filename}")
    if os.path.exists(file_path):
        new_filename = expected_filename.replace("_temp_folder.txt", ".txt")
        new_file_path = os.path.join(outdir, new_filename)

        try:
            os.rename(file_path, new_file_path)
        except Exception as e:
            print(f"Error renaming file : {e}")
    else:
        print(f"The {file_path} file does not exist.")


def move_file_to_outdir_folder(temporary_folder: str, outdir_folder: str) -> None:
    """
    Move files from temporary directory to final output directory.

    Used in overwrite mode to replace existing output with new results.
    Clears the existing output directory, renames the results file,
    and transfers all files from the temporary directory.

    Parameters
    ----------
    temporary_folder : str
        Path to the temporary processing directory.
    outdir_folder : str
        Path to the final output directory.
    """
    for fichier in os.listdir(outdir_folder):
        chemin_fichier = os.path.join(outdir_folder, fichier)
        try:
            if os.path.isfile(chemin_fichier):
                os.unlink(chemin_fichier)
            elif os.path.isdir(chemin_fichier):
                shutil.rmtree(chemin_fichier)
        except Exception as e:
            print(f"Error when deleting {chemin_fichier}: {e}")

    rename_temp_folder_file(outdir_folder, temporary_folder)

    for fichier in os.listdir(temporary_folder):
        source_path = os.path.join(temporary_folder, fichier)
        destination_path = os.path.join(outdir_folder, fichier)
        try:
            if os.path.isfile(source_path):
                shutil.move(source_path, destination_path)
            elif os.path.isdir(source_path):
                shutil.move(source_path, destination_path)
        except Exception as e:
            print(f"Error when transferring {source_path} to {destination_path}:: {e}")

    try:
        shutil.rmtree(temporary_folder)
    except Exception as e:
        print(f"Error when deleting the temporary folder {temporary_folder}: {e}")


def create_common_parser():
    """
    Create a parser with arguments common to multiple subcommands.

    Returns
    -------
    argparse.ArgumentParser
        Parser with common arguments for assembly analysis commands.
    """
    common = argparse.ArgumentParser(add_help=False)

    common.add_argument('-a', '--assemblies', nargs='+', type=str, required=True,
                        help='FASTA file(s) for assemblies.')

    common.add_argument('-o', '--outdir', type=str,
                        default="results_" + datetime.datetime.today().strftime("%Y-%m-%d_%I-%M-%S_%p"),
                        help='Folder for detailed output (default: results_YYYY-MM-DD_II-MM-SS_PP)')

    common.add_argument('--min_identity', type=float, default=80.0,
                        help='Minimum alignment identity for main results (default: 80)')

    common.add_argument('--min_coverage', type=float, default=50.0,
                        help='Minimum alignment coverage for main results (default: 50)')

    common.add_argument('--threads', type=int, default=4,
                        help='The number of threads to use for processing. (default: 4)')

    common.add_argument('--overwrite', action='store_true',
                        help='Allows the output directory to be overwritten if it already exists')

    common.add_argument('--format', type=str, default='tsv',
                        choices=['tsv', 'json', 'html', 'all'],
                        help='Output format: tsv, json, html, or all (default: tsv)')

    common.add_argument('--no-qc', action='store_true',
                        help='Skip assembly quality control metrics calculation')

    # Verbose/quiet mode (mutually exclusive)
    verbosity = common.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true',
                           help='Enable verbose output with timestamps and debug information')
    verbosity.add_argument('-q', '--quiet', action='store_true',
                           help='Suppress all output except errors')

    return common


def parse_arguments():
    """
    Parse and validate command-line arguments for diphtOscan.

    Supports subcommand-based architecture with backward compatibility.
    Subcommands: all, species, mlst, amr, tox, qc, update

    Returns
    -------
    argparse.Namespace
        Parsed arguments with the following key attributes:
        - command: Subcommand name ('all', 'species', 'mlst', etc.)
        - assemblies: List of input FASTA file paths
        - outdir: Output directory path
        - format: Output format ('tsv', 'json', 'both')
        - And various boolean flags for analyses

    Raises
    ------
    SystemExit
        If no arguments are provided or if required arguments are missing.
    """
    parser = argparse.ArgumentParser(
        prog='diphtoscan',
        description='diphtOscan is a tool to screen genome assemblies '
                    'of the diphtheriae species complex (CdSC)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  diphtoscan all -a genome.fasta -st -t -res_vir    Run full analysis
  diphtoscan species -a genome.fasta                Species identification only
  diphtoscan mlst -a genome.fasta                   MLST typing only
  diphtoscan amr -a genome.fasta                    AMR/virulence screening
  diphtoscan tox -a genome.fasta                    Tox allele detection
  diphtoscan qc -a genome.fasta                     Assembly QC metrics only
  diphtoscan update                                 Update databases
        '''
    )

    parser.add_argument('--version', action='version', version='diphtOscan v' + __version__,
                        help="Show program's version number and exit")

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Common parser for shared arguments
    common = create_common_parser()

    # --- 'all' subcommand (full analysis) ---
    sp_all = subparsers.add_parser('all', parents=[common],
                                   help='Run complete analysis pipeline')

    sp_all.add_argument('-st', '--mlst', action='store_true',
                        help='Turn on species CdSC and MLST sequence type (default: no)')

    sp_all.add_argument('-t', '--tox', action='store_true',
                        help='Turn on tox allele (default: no)')

    sp_all.add_argument('-res_vir', '--resistance_virulence', action='store_true',
                        help='Turn on resistance and main virulence genes screening (default: no)')

    sp_all.add_argument('-plus', '--extend_genotyping', action='store_true',
                        help='Turn on all virulence genes screening (default: no)')

    sp_all.add_argument('-integron', '--integron', action='store_true',
                        help='Screening the integron (default: no)')

    sp_all.add_argument('-tree', '--tree', action='store_true',
                        help='Generates a phylogenetic tree from JolyTree')

    # --- 'species' subcommand ---
    sp_species = subparsers.add_parser('species', parents=[common],
                                       help='Species identification only')

    # --- 'mlst' subcommand ---
    sp_mlst = subparsers.add_parser('mlst', parents=[common],
                                    help='MLST typing only')

    # --- 'amr' subcommand ---
    sp_amr = subparsers.add_parser('amr', parents=[common],
                                   help='AMR and virulence gene screening')

    sp_amr.add_argument('-plus', '--extend_genotyping', action='store_true',
                        help='Turn on all virulence genes screening (default: no)')

    sp_amr.add_argument('-integron', '--integron', action='store_true',
                        help='Screening the integron (default: no)')

    # --- 'tox' subcommand ---
    sp_tox = subparsers.add_parser('tox', parents=[common],
                                   help='Tox allele detection only')

    # --- 'qc' subcommand ---
    sp_qc = subparsers.add_parser('qc', parents=[common],
                                  help='Assembly QC metrics only')

    # --- 'update' subcommand ---
    sp_update = subparsers.add_parser('update',
                                      help='Update MLST, Tox Allele & AMR databases')
    # Add verbosity options to update command
    update_verbosity = sp_update.add_mutually_exclusive_group()
    update_verbosity.add_argument('-v', '--verbose', action='store_true',
                                  help='Enable verbose output')
    update_verbosity.add_argument('-q', '--quiet', action='store_true',
                                  help='Suppress all output except errors')

    # Handle backward compatibility: if first arg is not a subcommand
    # and -a or --assemblies is present, assume 'all' command
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    # Check for legacy invocation (no subcommand)
    known_commands = ['all', 'species', 'mlst', 'amr', 'tox', 'qc', 'update']
    # Flags that argparse should handle directly (not legacy mode)
    argparse_flags = ['-h', '--help', '--version']
    if len(sys.argv) > 1 and sys.argv[1] not in known_commands and not sys.argv[1].startswith('-'):
        # First arg doesn't look like a command or flag, show help
        parser.print_help(file=sys.stderr)
        sys.exit(1)
    elif len(sys.argv) > 1 and sys.argv[1].startswith('-') and sys.argv[1] not in argparse_flags:
        # Legacy mode: flags without subcommand (but not --help/--version)
        if '-u' in sys.argv or '--update' in sys.argv:
            sys.argv.insert(1, 'update')
        elif '-a' in sys.argv or '--assemblies' in sys.argv:
            sys.argv.insert(1, 'all')
        else:
            parser.print_help(file=sys.stderr)
            sys.exit(1)

    args = parser.parse_args()
    args.extract = False

    args.path = os.path.dirname(os.path.abspath(__file__))

    # Set default values for flags not present in all subcommands
    if not hasattr(args, 'mlst'):
        args.mlst = args.command == 'mlst'
    if not hasattr(args, 'tox'):
        args.tox = args.command == 'tox'
    if not hasattr(args, 'resistance_virulence'):
        args.resistance_virulence = args.command == 'amr'
    if not hasattr(args, 'extend_genotyping'):
        args.extend_genotyping = False
    if not hasattr(args, 'integron'):
        args.integron = False
    if not hasattr(args, 'tree'):
        args.tree = False
    if not hasattr(args, 'assemblies'):
        args.assemblies = None
    if not hasattr(args, 'no_qc'):
        args.no_qc = False
    if not hasattr(args, 'format'):
        args.format = 'tsv'
    if not hasattr(args, 'verbose'):
        args.verbose = False
    if not hasattr(args, 'quiet'):
        args.quiet = False

    if args.command != 'update':
        args = test_required_dependency(args)

    return args


def write_json_output(results: pd.DataFrame, outdir: str, version: str) -> str:
    """
    Write analysis results to a JSON file.

    Parameters
    ----------
    results : pd.DataFrame
        Results DataFrame with strains as rows.
    outdir : str
        Output directory path.
    version : str
        diphtOscan version string.

    Returns
    -------
    str
        Path to the written JSON file.
    """
    # Convert DataFrame to dictionary
    samples_dict = results.to_dict(orient='index')

    # Build summary statistics
    summary = {
        'total_samples': len(results),
    }

    if 'species' in results.columns:
        summary['species_counts'] = results['species'].value_counts().to_dict()

    if 'ST' in results.columns:
        st_counts = results['ST'].value_counts().to_dict()
        summary['st_counts'] = st_counts

    if 'biovar' in results.columns:
        summary['biovar_counts'] = results['biovar'].value_counts().to_dict()

    if 'tox_allele' in results.columns:
        tox_positive = len(results[results['tox_allele'] != '-'])
        summary['tox_positive'] = tox_positive
        summary['tox_negative'] = len(results) - tox_positive

    # Build JSON structure
    json_output = {
        'version': version,
        'date': datetime.datetime.now().isoformat(),
        'samples': samples_dict,
        'summary': summary
    }

    # Write JSON file
    json_path = os.path.join(outdir, f"{os.path.basename(outdir)}.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(json_output, f, indent=2, default=str)

    return json_path


def run_species_analysis(genome: str, args) -> dict:
    """Run species identification for a single genome."""
    return get_species_results(genome, args.path + '/data/species', str(args.threads))


def run_mlst_analysis(genome: str, dict_genome: dict, MLST_db: tuple, args) -> dict:
    """Run MLST typing for a single genome."""
    cd_complex = is_cd_complex(dict_genome)
    return get_chromosome_mlst_results(MLST_db, genome, cd_complex, args)


def run_tox_analysis(genome: str, TOX_db: tuple, args) -> dict:
    """Run tox allele detection for a single genome."""
    return get_tox_results(TOX_db, genome, args)


def run_qc_analysis(genome: str) -> dict:
    """Run assembly QC analysis for a single genome."""
    try:
        stats = calculate_assembly_stats(genome)
        return format_qc_results(stats)
    except Exception as e:
        log_warning(f"Could not calculate QC metrics for {genome}: {e}")
        return {}


def run_amr_analysis(genome: str, strain: str, args, resistance_db: str) -> pd.DataFrame:
    """
    Run AMR/virulence screening for a single genome.

    Returns
    -------
    pd.DataFrame
        AMRFinderPlus results for this genome.
    """
    min_identity = "-1"  # Default amrfinder
    subprocess.run([
        'amrfinder',
        '--nucleotide', genome,
        '--name', strain,
        '--nucleotide_output', args.outdir + "/" + strain + ".prot.fa",
        '--output', args.outdir + "/" + strain + ".blast.out",
        '--ident_min', min_identity,
        '--coverage_min', str(args.min_coverage/100),
        '--organism', 'Corynebacterium_diphtheriae',
        '--database', resistance_db,
        '--threads', str(args.threads),
        '--translation_table', '11',
        '--plus',
        '--quiet'
    ], check=False)

    if is_non_zero_file(args.outdir + '/' + strain + ".prot.fa"):
        data = pd.read_csv(args.outdir + '/' + strain + ".blast.out", sep="\t", dtype='str')
        data['File'] = genome
        return data
    else:
        prot_fa = os.path.join(args.outdir, strain + ".prot.fa")
        blast_out = os.path.join(args.outdir, strain + ".blast.out")
        if os.path.exists(prot_fa):
            os.remove(prot_fa)
        if os.path.exists(blast_out):
            os.remove(blast_out)
        return pd.DataFrame()


def run_integron_analysis(genome: str, strain: str, args) -> dict:
    """Run integron detection for a single genome."""
    subprocess.run([
        'integron_finder',
        '--cpu', str(args.threads),
        '--outdir', args.outdir + "/",
        '--gbk',
        '--func-annot',
        '--mute',
        genome
    ], check=False)

    # Remove empty directories from integron finder results
    for dir_path in glob.glob(args.outdir + "/Results_Integron_Finder_*/"):
        for root, dirs, files_list in os.walk(dir_path, topdown=False):
            for d in dirs:
                dir_to_check = os.path.join(root, d)
                if os.path.isdir(dir_to_check) and not os.listdir(dir_to_check):
                    os.rmdir(dir_to_check)

    summary_path = args.outdir + "/Results_Integron_Finder_" + strain + "/" + strain + ".summary"
    try:
        files = pd.read_csv(summary_path, sep="\t", index_col=0, skiprows=2)
        return files[['CALIN', 'complete', 'In0']].sum().to_dict()
    except Exception:
        return {'CALIN': 0, 'complete': 0, 'In0': 0}


def main():
    """
    Main entry point for diphtOscan.

    Orchestrates the complete analysis pipeline:
    1. Parse command-line arguments and validate dependencies
    2. Update databases if requested (update subcommand)
    3. For each input assembly:
       - Calculate QC metrics (unless --no-qc)
       - Identify species using Mash distance
       - Perform MLST typing (if mlst subcommand or --mlst)
       - Detect tox alleles (if tox subcommand or --tox)
       - Screen for AMR and virulence genes (if amr subcommand or --resistance_virulence)
       - Detect integrons (if --integron)
       - Determine biovar (when applicable)
    4. Generate combined results table
    5. Create iTOL visualization templates
    6. Build phylogenetic tree (if --tree and >= 4 assemblies)
    7. Export results in requested format (TSV, JSON, or both)

    The function writes results to the specified output directory and exits
    with status 0 on success.

    Raises
    ------
    SystemExit
        On missing dependencies, invalid input, or file system errors.
    """
    args = parse_arguments()

    # Initialize logging based on verbosity flags
    logger = setup_logging(verbose=args.verbose, quiet=args.quiet)
    log_info(f"diphtOscan v{__version__}")

    MLST_db = (get_chromosome_mlst_header(), args.path + '/data/mlst/pubmlst_diphtheria_seqdef_scheme_3.fas', args.path + '/data/mlst/st_profiles.txt')
    TOX_db = (get_tox_header(), args.path + '/data/tox/pubmlst_diphtheria_seqdef_scheme_4.fas', args.path + '/data/tox/tox_profiles.txt')

    # Handle update command
    if args.command == 'update':
        update_database(args, MLST_db, TOX_db)
        sys.exit(0)

    update_database(args, MLST_db, TOX_db)

    if args.assemblies is None:
        sys.exit(0)

    resistance_db = find_resistance_db(args)

    if args.overwrite:
        args, final_output_path = redefine_output_file(args)

    try:
        os.makedirs(args.outdir)
        log_debug(f"Directory '{args.outdir}' created successfully")
    except OSError:
        log_error(f"Directory '{args.outdir}' can not be created")
        sys.exit(0)

    dict_results = {}
    data_resistance = pd.DataFrame()

    # Set up progress bar for multiple assemblies
    assemblies_iter = args.assemblies
    if not args.quiet and len(args.assemblies) > 1:
        assemblies_iter = tqdm(
            args.assemblies,
            desc="Processing",
            unit="file",
            file=sys.stderr
        )

    for genome in assemblies_iter:
        log_debug(f"Processing file: {genome}")
        basename = os.path.basename(genome)
        strain = os.path.splitext(basename)[0]

        dict_genome = {}

        # QC analysis (unless disabled)
        if not args.no_qc:
            qc_results = run_qc_analysis(genome)
            dict_genome.update(qc_results)

        # Species identification (always run for context)
        species_results = run_species_analysis(genome, args)
        dict_genome.update(species_results)

        # MLST typing
        if args.mlst or args.command in ['mlst', 'all']:
            if args.command == 'all' and not args.mlst:
                pass  # Skip if 'all' but --mlst not specified
            else:
                mlst_results = run_mlst_analysis(genome, dict_genome, MLST_db, args)
                dict_genome.update(mlst_results)

        # Tox allele detection
        if args.tox or args.command in ['tox', 'all']:
            if args.command == 'all' and not args.tox:
                pass  # Skip if 'all' but --tox not specified
            else:
                tox_results = run_tox_analysis(genome, TOX_db, args)
                dict_genome.update(tox_results)

        # AMR/Virulence screening
        if args.resistance_virulence or args.command == 'amr':
            amr_data = run_amr_analysis(genome, strain, args, resistance_db)
            if not amr_data.empty:
                data_resistance = pd.concat([data_resistance, amr_data], axis=0, ignore_index=True)
                dict_genome.update({"GENOMIC_CONTEXT": get_genomic_context(args.outdir, amr_data)})

        # Integron detection
        if args.integron:
            integron_results = run_integron_analysis(genome, strain, args)
            dict_genome.update(integron_results)

        # Biovar determination (after all gene detection is complete)
        species = dict_genome.get('species', '')
        biovar = determine_biovar(dict_genome, species)
        dict_genome['biovar'] = biovar

        dict_results[strain] = dict_genome

    table_results = pd.DataFrame(dict_results)
    table_results = table_results.T

    if len(data_resistance.index) != 0:
        table_resistance = armfinder_to_table(data_resistance)
        for family in table_resistance.columns:
            table_resistance[family] = table_resistance[family].apply(lambda x: ";".join(sorted(x.split(';'))))

        if not args.extend_genotyping:
            header = [x for x in table_resistance.columns if x not in delete_virulence_extended()]
            table_resistance = table_resistance[sorted(header)]

        table_resistance = table_resistance.replace('', '-')
        results = pd.concat([table_results, table_resistance], axis=1, join='outer')
    else:
        results = table_results

    results = results.infer_objects().fillna("-")

    # Generate iTOL templates
    spuA(results, args)
    narG(results, args)
    toxin(results, args)
    amr_families(results, args)

    # Calculate summary statistics
    summary = calculate_summary_statistics(results)

    # Print summary to console (unless quiet mode)
    if not args.quiet:
        summary_text = format_summary_for_console(summary)
        log_info(summary_text)

    # Write output files based on format
    output_basename = args.outdir.split("/")[-1]

    if args.format in ('tsv', 'all'):
        results.to_csv(args.outdir + "/" + output_basename + ".txt", sep='\t')
        log_info(f"Results written to: {args.outdir}/{output_basename}.txt")

    if args.format in ('json', 'all'):
        json_path = write_json_output(results, args.outdir, __version__)
        log_info(f"JSON results written to: {json_path}")

    if args.format in ('html', 'all'):
        html_path = generate_html_report(results, summary, args.outdir, __version__)
        log_info(f"HTML report: {html_path}")

    # Generate phylogenetic tree
    if hasattr(args, 'tree') and args.tree and len(args.assemblies) >= 4:
        generate_jolytree(args)

    if args.overwrite:
        move_file_to_outdir_folder(temporary_folder=args.outdir,
                                   outdir_folder=final_output_path)


if __name__ == "__main__":
    main()
