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
    python diphtOscan.py argument1 argument2

    argument1: un entier signifiant un truc
    argument2: une chaîne de caractères décrivant un bidule
"""

__authors__ = ("Melanie HENNART; Martin RETHORET-PASTY")
__contact__ = ("martin.rethoret-pasty@pasteur.fr")
__version__ = "1.7.0" 
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
import pandas as pd
import datetime
import os.path
import argparse
import subprocess
import shutil

sys.path.append('/module')

from typing import List
from module.species import get_species_results, is_cd_complex
from module.template_iTOL import spuA, narG, toxin, amr_families
from module.updating_database import update_database
from module.jolytree_generation import generate_jolytree

from module.utils import (
    get_chromosome_mlst_results, 
    get_tox_results, 
    get_chromosome_mlst_header, 
    get_tox_header,
    delete_virulence_extended,
    is_non_zero_file,
    armfinder_to_table,
    get_genomic_context,
    find_resistance_db
    )

def test_unique_dependency(name:str):
    return shutil.which(name) is not None


def test_multiple_dependencies(dependencies:List[str]):
    for dependency in dependencies:
        presence = test_unique_dependency(dependency)
        if presence is not True:
            print(f'/!\\ Warning /!\\ : {dependency} missing in path!')
            sys.exit(-1)


def test_required_dependency(args):
    diphtoscan_dependencies = ["mash",'amrfinder','hmmsearch', 'makeblastdb','blastn', 'blastp'] 
    joly_tree_path = "diphtoscan/script/JolyTree/JolyTree.sh"
    joly_tree_dependencies = ["gawk",'fastme','REQ']
    integron_fender_dependencies = ['hmmsearch', 'cmsearch', 'prodigal']
    
    if args.assemblies == None: #TODO :Ensure that dependencies are not required to update the database
        return args

    subprocess.run(["echo", "Dependency testing"])
    test_multiple_dependencies(diphtoscan_dependencies)
    
    if args.integron: 
        rc = test_unique_dependency("Integron_finder")
        test_multiple_dependencies(integron_fender_dependencies)
        if rc == 0:
            args.integron = True
        else:
            print('/!\\ Warning /!\\ : Integron_finder missing in path! Integron analysis not carried out.')
            args.integron = False

    if args.tree:
        if os.path.isfile(joly_tree_path):
            test_multiple_dependencies(joly_tree_dependencies)
            args.tree = True
        else:
            print('/!\\ Warning /!\\ : JolyTree.sh missing in /diphtoscan/script/JolyTree/ ! Joly_tree representation not carried out.')
            args.tree = False
    print('\n')
    return args


def redefine_output_file(args):

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


def rename_temp_folder_file(directory):
    path, file_name = os.path.split(directory)
    expected_filename = f"{file_name}.txt"
    file_path = os.path.join(directory, expected_filename)

    if os.path.exists(file_path):
        new_filename = expected_filename.replace("_temp_folder.txt", ".txt")
        new_file_path = os.path.join(args.outdir, new_filename)
        
        try:
            os.rename(file_path, new_file_path)
        except Exception as e:
            print(f"Error renaming file : {e}")
    else:
        print(f"The {file_path} file does not exist.")


def move_file_to_outdir_folder(temporary_folder, outdir_folder):
    
    for fichier in os.listdir(outdir_folder):
        chemin_fichier = os.path.join(outdir_folder, fichier)
        try:
            if os.path.isfile(chemin_fichier):
                os.unlink(chemin_fichier)
            elif os.path.isdir(chemin_fichier):
                shutil.rmtree(chemin_fichier)
        except Exception as e:
            print(f"Error when deleting {chemin_fichier}: {e}")

    rename_temp_folder_file(temporary_folder)

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


def parse_arguments():
    parser = argparse.ArgumentParser(description='diphtOscan is a tool to screen genome assemblies '
                                                 'of the diphtheriae species complex (CdSC)',
                                     add_help=False)

    updating_args = parser.add_argument_group('Updating option')

    updating_args.add_argument('-u', '--update', action='store_true',
                                help='Update database MLST, Tox Allele & AMR (default: no).'
                                'The database update can be executed on its own without the -a option.')
    
    required_args = parser.add_argument_group('Required option')
    required_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=('-u' not in sys.argv and '--update' not in sys.argv),
                               help='FASTA file(s) for assemblies. ') #-a is required only if -u is not present. It allows the user to update the database easily

    screening_args = parser.add_argument_group('Screening options')
                             
    screening_args.add_argument('-st', '--mlst', action='store_true',
                                help='Turn on species Corynebacterium diphtheriae species complex (CdSC)'
                                     ' and MLST sequence type (default: no)')

    screening_args.add_argument('-t', '--tox', action='store_true',
                                help='Turn on tox allele (default: no)')

    screening_args.add_argument('-res_vir', '--resistance_virulence', action='store_true',
                                help='Turn on resistance and main virulence genes screening (default: no resistance '
                                     'and virulence gene screening)')

    screening_args.add_argument('-plus', '--extend_genotyping', action='store_true',
                                help='Turn on all virulence genes screening (default: no all virulence '
                                     'gene screening)')

    screening_args.add_argument('-integron', '--integron', action='store_true',
                                help='Screening the intregon(default: no)')
                                     
    output_args = parser.add_argument_group('Output options')

     
    output_args.add_argument('-o', '--outdir', type=str, default="results_"+ datetime.datetime.today().strftime("%Y-%m-%d_%I-%M-%S_%p"),
                             help='Folder for detailed output (default: results_YYYY-MM-DD_II-MM-SS_PP)')

    setting_args = parser.add_argument_group('Settings')
    
    setting_args.add_argument('--min_identity', type=float, default=80.0,
                              help='Minimum alignment identity for main results (default: 80)')
    
    setting_args.add_argument('--min_coverage', type=float, default=50.0,
                              help='Minimum alignment coverage for main results (default: 50)')
    
    setting_args.add_argument('--threads', type=int, default=4,
                              help='The number of threads to use for processing. (default: 4)')
    
    setting_args.add_argument('--overwrite', action='store_true',
                              help='Allows the output directory to be overwritten if it already exists')
    
    tree_args = parser.add_argument_group('Phylogenetic tree')
    tree_args.add_argument('-tree', '--tree', action='store_true',
                           help='Generates a phylogenetic tree from JolyTree')
    
    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='diphtOscan v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the entire help (argparse default is to just give an error
    # like '-a is required').
    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.extract = False    
    args.path = os.path.dirname(os.path.abspath(__file__))

    args = test_required_dependency(args)
    return args 


if __name__ == "__main__":
      
    args = parse_arguments()
    get_path = os.getcwd()

    MLST_db = (get_chromosome_mlst_header(), args.path + '/data/mlst/pubmlst_diphtheria_seqdef_scheme_3.fas', args.path + '/data/mlst/st_profiles.txt') 
    TOX_db = (get_tox_header(), args.path + '/data/tox/pubmlst_diphtheria_seqdef_scheme_4.fas', args.path + '/data/tox/tox_profiles.txt') 

    update_database(args,MLST_db,TOX_db)
    
    if args.assemblies == None:
        sys.exit(0)

    resistance_db = find_resistance_db(args) 
    
    if args.overwrite :
        args, final_output_path = redefine_output_file(args)

    try:
        os.makedirs(args.outdir)
        print("Directory '%s' created successfully \n" %args.outdir)
    except OSError :
        print("Directory '%s' can not be created \n"  %args.outdir)        
        sys.exit(0)
	
    dict_results = {}
    data_resistance = pd.DataFrame()
    for genome in args.assemblies :

        basename = os.path.basename(genome)
        strain = os.path.splitext(basename)[0]

        dict_genome =  get_species_results(genome, args.path + '/data/species', str(args.threads)) 
        if args.mlst : 
            cd_complex = is_cd_complex(dict_genome)
            dict_genome.update(get_chromosome_mlst_results(MLST_db, genome, cd_complex, args))
        
        if args.tox :
            dict_genome.update(get_tox_results(TOX_db, genome, args))
            
        if args.resistance_virulence :   
            min_identity = "-1" # Defaut amrfinder
            os.system('amrfinder --nucleotide ' + genome +
                      ' --name '+strain+
                      ' --nucleotide_output ' + args.outdir + "/" + strain + ".prot.fa" +
                      ' --output '+ args.outdir + "/" + strain + ".blast.out" +
                      ' --ident_min '+ min_identity +
                      ' --coverage_min ' + str(args.min_coverage/100) +
                      ' --organism Corynebacterium_diphtheriae' +
                      ' --database ' + resistance_db +
                      ' --threads ' + str(args.threads)+
                      #' --blast_bin /opt/gensoft/exe/blast+/2.12.0/bin/' +
                      ' --translation_table 11 --plus --quiet ')
            if is_non_zero_file(args.outdir +'/' +strain + ".prot.fa"):
                data = pd.read_csv(args.outdir +'/' + strain + ".blast.out",sep="\t", dtype='str')
                data['File'] = genome
                data_resistance = pd.concat([data_resistance, data], axis = 0, ignore_index=True)
                dict_genome.update({"GENOMIC_CONTEXT" : get_genomic_context (args.outdir, data)})
            else :
                os.system('rm '+ args.outdir +'/' + strain + ".prot.fa")
                os.system('rm '+ args.outdir +'/' + strain + ".blast.out")
                
        if args.integron :
          os.system('integron_finder --cpu ' + str(args.threads)+
                    ' --outdir '+ args.outdir + "/" +
                    ' --gbk --func-annot --mute '+ genome)   
          os.system('find '+ args.outdir + "/Results_Integron_Finder_*/ " + '-empty -type d -delete')

          files = pd.read_csv(args.outdir + "/Results_Integron_Finder_"+strain + "/" + strain+".summary",sep="\t", index_col=0, skiprows = 2)
          dict_genome.update(files[['CALIN','complete','In0']].sum().to_dict())
          
        dict_results[strain] = dict_genome
     
    table_results = pd.DataFrame(dict_results)
    table_results = table_results.T
    
    if len(data_resistance.index) != 0 :
        table_resistance = armfinder_to_table(data_resistance)
        for family in table_resistance.columns:
            table_resistance[family] = table_resistance[family].apply(lambda x : ";".join(sorted(x.split(';'))))

        if not args.extend_genotyping : 
            header = [x for x in table_resistance.columns if x not in delete_virulence_extended()] 
            table_resistance = table_resistance[sorted(header)]

        table_resistance = table_resistance.replace('','-')    
        results = pd.concat([table_results, table_resistance], axis=1, join='outer')
    else : 
        results = table_results
        
    results = results.infer_objects().fillna("-")
        
    spuA(results, args)
    narG(results, args)
    toxin(results, args)
    amr_families(results, args)
    
    results.to_csv(args.outdir+"/"+args.outdir.split("/")[-1]+".txt", sep='\t')
    
    if args.tree and len(args.assemblies) >= 4 :
        generate_jolytree(args)
  
    if args.overwrite :
        move_file_to_outdir_folder(temporary_folder = args.outdir,
                                   outdir_folder = final_output_path
                                   )
