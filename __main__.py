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
     - Species (e.g. C. diphtheriae, C. belfantii, C. rouxii, C. ulcerans 
                and C. pseudotuberculosis)
     - MLST sequence type
     - Virulence loci 
     - Antimicrobial resistance: acquired genes, SNPs
     - Biovar prediction
     - Detection of tox gene 

Usage:
======
    python diphtOscan.py argument1 argument2

    argument1: un entier signifiant un truc
    argument2: une chaîne de caractères décrivant un bidule
"""

__authors__ = ("Melanie HENNART")
__contact__ = ("melanie.hennart@pasteur.fr")
__version__ = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2022/10/11"

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
import glob
import os.path
import networkx as nwx
import itertools
import argparse

sys.path.append('/module')

from module.download_alleles_st import create_db, download_profiles_st, download_profiles_tox
from module.species import get_species_results, is_cd_complex
from module.mlstBLAST import mlst_blast

 

         
    
def get_chromosome_mlst_results(infoMLST, contigs, cd_complex, args):
    chromosome_mlst_header = infoMLST[0]
    if cd_complex:
        seqs = infoMLST[1]
        database = infoMLST[2]
        chr_st, chr_st_detail, _, _ = \
             mlst_blast(seqs, database, 'no', [contigs], min_cov=args.min_coverage,
                       min_ident=args.min_identity, max_missing=3, allow_multiple=False)
        if chr_st != '0':
            chr_st = 'ST' + chr_st
        
        assert len(chromosome_mlst_header) == len(chr_st_detail)
        results = {'ST': chr_st}

    else:
        results = {'ST': "NA"}
        chr_st_detail = ['-'] * len(chromosome_mlst_header)

    results.update(dict(zip(infoMLST[0], chr_st_detail)))
    return results

def get_tox_results(infoTOX, contigs, args):
    tox_header = infoTOX[0]
    seqs = infoTOX[1]
    database = infoTOX[2]
    chr_st, chr_st_detail, _, _ = \
         mlst_blast(seqs, database, 'no', [contigs], min_cov=args.min_coverage,
                   min_ident=args.min_identity, max_missing=3, allow_multiple=False)
    if chr_st != '0':
        chr_st = 'TOX' + chr_st
    
    assert len(tox_header) == len(chr_st_detail)
    results = dict(zip(infoTOX[0], chr_st_detail))
    #results = {'TOX': chr_st}
    
    #results.update(dict(zip(infoTOX[0], chr_st_detail)))
    return results

def get_chromosome_mlst_header():
    return ['atpA', 'dnaE', 'dnaK', 'fusA', 'leuA', 'odhA', 'rpoB']

def get_tox_header():
    return ['tox']

def get_virulence():
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'spuA', 'narG',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN', 'iusABCDE',
            'chtAB','htaA-hmuTUV-htaBC', 'cdtQP-sidBA-ddpABCD']
    
def get_virulence_extended(): 
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'SpuA-CLUSTER', 'narIJHK',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','cdtQP-sidBA-ddpABCD']
def delete_virulence_extended():
    return [ 'SpuA-CLUSTER','narIJHK','SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','cdtQP-sidBA-ddpABCD']
       
       
def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def armfinder_to_table(data_resistance):
    dico_Method = {'ALLELEX' : "",
                   'EXACTX' :  "",
                   'POINTX' : "!",
                   'BLASTX' : "*",
                   'PARTIALX' : "?",
                   'PARTIAL_CONTIG_ENDX' : "?$",
                   'INTERNAL_STOP' :  "#"}
    data_resistance['Class'] = data_resistance['Class'].fillna ('NoClass')
    Class = data_resistance['Class'].value_counts().keys()
    Strains = data_resistance['Name'].value_counts().keys()
    table = pd.DataFrame('',index=Strains, columns=Class)
    for res in data_resistance.index : 
        gene = data_resistance['Gene symbol'][res] + dico_Method[data_resistance['Method'][res]]      
        if 'tox' in data_resistance['Gene symbol'][res] :    
            if float(data_resistance['% Coverage of reference sequence'][res]) != 100.00 : 
                if (data_resistance['Method'][res] == 'BLASTX') :
                    gene = data_resistance['Gene symbol'][res] + "-NTTB?-"+str(round(100-float(data_resistance['% Coverage of reference sequence'][res])))+"%"              
                else :
                    gene = data_resistance['Gene symbol'][res] + "-NTTB" + dico_Method[data_resistance['Method'][res]]             
        if (data_resistance['Method'][res] == 'PARTIALX') or \
           (data_resistance['Method'][res] == 'PARTIAL_CONTIG_ENDX') or \
           (data_resistance['Method'][res] == 'INTERNAL_STOP') : 
               gene += "-"+str(round(100-float(data_resistance['% Coverage of reference sequence'][res])))+"%"
        
        strain = data_resistance['Name'][res]
        family = data_resistance['Class'][res]
        
        if table[family][strain] != '' :
               table[family][strain] += ";"
        table[family][strain] += gene
    return table

def get_genomic_context (data) :
    d = []
    data_AMR = data[~data['Class'].isin( list(set(get_virulence_extended())| set(get_virulence())))]
    for contigs in data_AMR['Contig id'].value_counts().keys() :            
        table_contigs  = data_AMR[data_AMR['Contig id'] == contigs]
        
        if len(table_contigs) == 1 : 
            d.append({table_contigs['Gene symbol'].value_counts().keys()[0]})
        else:  
            DG = nwx.Graph()
            for (x,y) in itertools.combinations( table_contigs.index, 2 ):         
                a = int(table_contigs['Start'][x]) + (int(table_contigs['Stop'][x]) - int(table_contigs['Start'][x]) / 2)
                b = int(table_contigs['Start'][y]) + (int(table_contigs['Stop'][y]) - int(table_contigs['Start'][y]) / 2)
                DG.add_node(table_contigs['Gene symbol'][x])
                DG.add_node(table_contigs['Gene symbol'][y])    
                if abs(b - a) <= 8000 :      
                    DG.add_edge(table_contigs['Gene symbol'][x], table_contigs['Gene symbol'][y], weight=abs(b - a))    
            d.extend([set(c) for c in nwx.connected_components(DG)])
            
    return " $$ ".join([';'.join(x) for x in d])

def find_resistance_db (args):
            files = [ name for name in glob.glob(args.path+'/data/resistance/*') if os.path.isdir(name) ]
            max_file = max(files, key = os.path.getctime)    
            return max_file 
        
def parse_arguments():
    parser = argparse.ArgumentParser(description='diphtOscan: a tool for characterising '
                                                 'virulence and resistance in Corynebacterium',
                                     add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True,
                               help='FASTA file(s) for assemblies')

    screening_args = parser.add_argument_group('Screening options')
    
    screening_args.add_argument('-u', '--update', action='store_true',
                                help='Update database MLST et AMR (default: no)')
    
    screening_args.add_argument('-t', '--taxonomy', action='store_true',
                                help='Turn on species Corynebacterium diphtheriae species complex (CdSC)'
                                     ' and MLST sequence type (default: no)')
    
    screening_args.add_argument('-res_vir', '--resistance_virulence', action='store_true',
                                help='Turn on resistance and virulence genes screening (default: no resistance '
                                     'and virulence gene screening)')
    screening_args.add_argument('-plus', '--extend_genotyping', action='store_true',
                                help='Turn on all virulence genes screening (default: no all virulence '
                                     'gene screening)')
    
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
    
    if args.taxonomy:
        args.mlst = True
    else :
        args.mlst = False
        
    args.extract = False    
    args.path = os.path.dirname(os.path.abspath(__file__))  
    
    return args 


 
if __name__ == "__main__":
      
    #=== Parameters
    args = parse_arguments()
    get_path = os.getcwd()
    #=== 
      
    MLST_db = (get_chromosome_mlst_header(), args.path + '/data/mlst/pubmlst_diphtheria_seqdef_scheme_3.fas', args.path + '/data/mlst/st_profiles.txt') 
    TOX_db = (get_tox_header(), args.path + '/data/tox/pubmlst_diphtheria_seqdef_scheme_4.fas', args.path + '/data/tox/tox_profiles.txt') 

    if args.update : 
        os.system("rm "+ MLST_db[1] + "* " + MLST_db[2] + "* ")                  
        print("Downloading MLST database")
        path_mlst_sequences, loci_mlst = create_db("pubmlst_diphtheria_seqdef", "3", args.path +"/data/mlst")
        download_profiles_st ("pubmlst_diphtheria_seqdef", "3", args.path +"/data/mlst", loci_mlst)
        print("   ... done \n")
        
        print("Downloading tox database")
        path_tox_sequences, loci_tox = create_db("pubmlst_diphtheria_seqdef", "4", args.path +"/data/tox")
        download_profiles_tox ("pubmlst_diphtheria_seqdef", "4", args.path +"/data/tox")
        print("   ... done \n")
        
        print ('bash ' + args.path + '/data/resistance/update_database_resistance.sh')
        os.system('bash ' + args.path + '/data/resistance/update_database_resistance.sh')
        print("   ... done \n\n\n")
    try:
        os.makedirs(args.outdir)
        print("Directory '%s' created successfully \n" %args.outdir)
    except OSError :
        print("Directory '%s' can not be created \n"  %args.outdir)        
        sys.exit(0)
    
    resistance_db = find_resistance_db(args) 
    prediction_db = args.path +"/data/virulence"
    
    dict_results = {}
    data_resistance = pd.DataFrame()
    for genome in args.assemblies :
        strain = genome.split('/')[-1].split('.')[0]
        print (strain)
        fasta =  get_path +'/'+genome
        dict_genome =  get_species_results(fasta, args.path + '/data/species', str(args.threads))   
        
        if args.mlst : 
            cd_complex = is_cd_complex(dict_genome)
            dict_genome.update(get_chromosome_mlst_results(MLST_db, fasta, cd_complex, args))
        
        dict_genome.update(get_tox_results(TOX_db, fasta, args))
            
        if args.resistance_virulence :   
            min_identity = "-1" # defaut amrfinder
            #min_coverage = "0.8" # defaut amrfinder                
            os.system('amrfinder --nucleotide ' + fasta +
                      ' --name '+strain+
                      ' --nucleotide_output ' + args.outdir + "/" + strain + ".prot.fa" +
                      ' --output '+ args.outdir + "/" + strain + ".blast.out" +
                      ' --ident_min '+ min_identity +
                      ' --coverage_min ' + str(args.min_coverage/100) +
                      ' --organism Corynebacterium_diphtheriae' +
                      ' --database ' + resistance_db +
                      ' --threads ' + str(args.threads)+
                      ' --blast_bin /opt/gensoft/exe/blast+/2.12.0/bin/' +
                      ' --translation_table 11 --plus --quiet ')
            if is_non_zero_file(args.outdir +'/' +strain + ".prot.fa"):
                data = pd.read_csv(args.outdir +'/' + strain + ".blast.out",sep="\t", dtype='str')
                data_resistance = pd.concat([data_resistance, data], axis = 0, ignore_index=True)
                if  args.extend_genotyping :
                    dict_genome.update({"GENOMIC_CONTEXT" : get_genomic_context (data)})
            else :
                os.system('rm '+ args.outdir +'/' + strain + ".prot.fa")
                os.system('rm '+ args.outdir +'/' + strain + ".blast.out")
            
            
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
        
    results = results.fillna("-")
    results.to_csv(args.outdir+"/"+args.outdir.split("/")[-1]+".txt", sep='\t')
    
    
    if args.tree and len(args.assemblies) >= 4 : 
        print ("\nGenerating a phylogenetic tree from JolyTree \n")
        os.makedirs(args.outdir+"/FolderJolyTree" )
        os.system("cp "+ " ".join(args.assemblies) + " " + args.outdir+"/FolderJolyTree/")   
        os.system('bash ' + args.path + '/script/JolyTree/JolyTree.sh -i '+ args.outdir+ "/FolderJolyTree -b "+ args.outdir+"/"+args.outdir + ".jolytree -t " + str(args.threads) + " > "
                  + args.outdir+"/"+args.outdir + ".jolytree.log")
        os.system("rm -R "+ args.outdir+"/FolderJolyTree/")
   


