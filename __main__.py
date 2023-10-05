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
__version__ = "1.2.0"
__copyright__ = "copyleft"
__date__ = "2023/04/12"

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
import argparse
import subprocess

sys.path.append('/module')

from module.species import get_species_results, is_cd_complex
from module.mlstBLAST import mlst_blast
from module.template_iTOL import writeTemplateStrip
from module.biovar_detection import spuA, narG, toxin
from module.updating_databse import update_database
from module.jolytree_generation import generate_jolytree

         
    
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

def get_chromosome_mlst_header()-> list:
    return ['atpA', 'dnaE', 'dnaK', 'fusA', 'leuA', 'odhA', 'rpoB']

def get_tox_header() -> list:
    return ['tox']

def get_virulence() -> list:
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'spuA', 'narG',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN', 'iusABCDE',
            'chtAB','htaA-hmuTUV-htaBC', 'cdtQP-sidBA-ddpABCD']
    
def get_virulence_extended() -> list: 
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'SpuA-CLUSTER', 'narIJHK',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']
def delete_virulence_extended() -> list:
    return [ 'SpuA-CLUSTER','narIJHK','SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']
       
       
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

def get_genomic_context (outdir, data) :
    d = []
    data_AMR = data[~data['Class'].isin( list(set(get_virulence_extended())| set(get_virulence())))]
    fi = open(outdir+'/distance_context.txt', 'a', encoding='utf-8')
    for contigs in data_AMR['Contig id'].value_counts().keys() :            
        table_contigs  = data_AMR[data_AMR['Contig id'] == contigs]
        
        if len(table_contigs) == 1 : 
            d.append(table_contigs['Gene symbol'].value_counts().keys()[0])
        else:  
            t = table_contigs['Gene symbol'].iloc[0]
            for i in range(0,len(table_contigs)-1):
                dis = int(table_contigs['Start'].iloc[i+1]) - int(table_contigs['Stop'].iloc[i])
                fi.write(table_contigs['Gene symbol'].iloc[i]+'\t'+table_contigs['Gene symbol'].iloc[i+1]+'\t'+str(abs(dis))+'\n')
                if abs(dis) <=  8000 :  
                    t +=  ";" + table_contigs['Gene symbol'].iloc[i+1]
                else :
                    t +=  " || " + table_contigs['Gene symbol'].iloc[i+1]                
            d.append(t)
    fi.close()        
    return " || ".join(d)

def find_resistance_db (args):
            files = [ name for name in glob.glob(args.path+'/data/resistance/*') if os.path.isdir(name) ]
            max_file = max(files, key = os.path.getctime)    
            return max_file 
        
def parse_arguments():
    parser = argparse.ArgumentParser(description='diphtOscan: a tool for characterising '
                                                 'virulence and resistance in Corynebacterium',
                                     add_help=False)

    screening_args = parser.add_argument_group('Screening options')
    
    screening_args.add_argument('-u', '--update', action='store_true',
                                help='Update database MLST et AMR (default: no)')
    
    screening_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=('-u' not in sys.argv and '--update' not in sys.argv),
                               help='FASTA file(s) for assemblies') #-a is required only if -u is not present. It allows the user to update the database easily
                                
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
    
    return args 


    if args.integron: 
        rc = subprocess.call(["command", "-v", "integron_finder"])
        if rc == 0:
            args.integron = True
        else:
            print('integron_finder missing in path!')
            args.integron = False
 
if __name__ == "__main__":
      
    args = parse_arguments()
    get_path = os.getcwd()

    MLST_db = (get_chromosome_mlst_header(), args.path + '/data/mlst/pubmlst_diphtheria_seqdef_scheme_3.fas', args.path + '/data/mlst/st_profiles.txt') 
    TOX_db = (get_tox_header(), args.path + '/data/tox/pubmlst_diphtheria_seqdef_scheme_4.fas', args.path + '/data/tox/tox_profiles.txt') 

    update_database(args)
    
    resistance_db = find_resistance_db(args) 
    prediction_db = args.path +"/data/virulence"
    
    dict_results = {}
    data_resistance = pd.DataFrame()
    for genome in args.assemblies :
        basename = os.path.basename(genome)
        strain = os.path.splitext(basename)[0]
        
        fasta =  get_path +'/'+genome
        dict_genome =  get_species_results(fasta, args.path + '/data/species', str(args.threads))   
        
        if args.mlst : 
            cd_complex = is_cd_complex(dict_genome)
            dict_genome.update(get_chromosome_mlst_results(MLST_db, fasta, cd_complex, args))
        
        if args.tox :
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
                      #' --blast_bin /opt/gensoft/exe/blast+/2.12.0/bin/' +
                      ' --translation_table 11 --plus --quiet ')
            if is_non_zero_file(args.outdir +'/' +strain + ".prot.fa"):
                data = pd.read_csv(args.outdir +'/' + strain + ".blast.out",sep="\t", dtype='str')
                data_resistance = pd.concat([data_resistance, data], axis = 0, ignore_index=True)
                dict_genome.update({"GENOMIC_CONTEXT" : get_genomic_context (args.outdir, data)})
            else :
                os.system('rm '+ args.outdir +'/' + strain + ".prot.fa")
                os.system('rm '+ args.outdir +'/' + strain + ".blast.out")
                
        if args.integron :
          os.system('integron_finder --cpu ' + str(args.threads)+
                    ' --outdir '+ args.outdir + "/" +
                    ' --gbk --func-annot --mute '+ fasta)   
          os.system('find '+ args.outdir + "/Results_Integron_Finder_*/ " + '-empty -type d -delete')
          
          files = pd.read_csv(args.outdir + "/Results_Integron_Finder_"+strain + "/" + strain+".summary",sep="\t", index_col=0, skiprows = 1)
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
        
    results = results.fillna("-")
    
    # spuA and narG 
    
    spuA(results, args)
    narG(results, args)
    toxin(results, args)
    
    # ARM
    list_familiesRes ={'AMINOGLYCOSIDE' : ['#a6cee3', '#1f78b4'],
                            'MACROLIDE' : ['#b2df8a', '#33a02c'],
                             'PHENICOL' : ['#fb9a99', '#e31a1c'],
                          'SULFONAMIDE' : ['#fdbf6f', '#ff7f00'],
                         'TETRACYCLINE' : ['#cab2d6', "#6a3d9a"],
                         'TRIMETHOPRIM' : ['#ffff99', '#b15928'],
                  'QUATERNARY AMMONIUM' : ["#e0eaf4", "#3c6498"],
                          'BETA-LACTAM' : ["#da74da", "#9e3c8b"],
                            'QUINOLONE' : ["#b0c665", "#6c7b38"],
                            'RIFAMYCIN' : ["#bd924f", "#926114"]}


    for family in list_familiesRes : 
        if family in results.columns:
            writeTemplateStrip (args.outdir, results, family, list_familiesRes)
    
    
    
    results.to_csv(args.outdir+"/"+args.outdir.split("/")[-1]+".txt", sep='\t')
    
    if args.tree and len(args.assemblies) >= 4 :
        generate_jolytree(args)
   


