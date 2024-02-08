import sys
import os

from module.download_alleles_st import create_db, download_profiles_st, download_profiles_tox

def update_database(arguments, mlst_database:tuple, tox_database:tuple):
    if arguments.update : 
        os.system("rm "+ mlst_database[1] + "* " + mlst_database[2] + "* ")                  
        print("Downloading MLST database")
        path_mlst_sequences, loci_mlst = create_db("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst")
        download_profiles_st ("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst", loci_mlst)
        print("   ... done \n")

        os.system("rm "+ tox_database[1] + "* " + tox_database[2] + "* ")                  
        print("Downloading tox database")
        path_tox_sequences, loci_tox = create_db("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        download_profiles_tox ("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        print("   ... done \n")
        
        os.system('bash ' + arguments.path + '/data/resistance/update_database_resistance.sh')
        print("   ... done \n\n\n")