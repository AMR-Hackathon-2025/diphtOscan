import sys
import os

def update_database(arguments):
    if arguments.update : 
        os.system("rm "+ MLST_db[1] + "* " + MLST_db[2] + "* ")                  
        print("Downloading MLST database")
        path_mlst_sequences, loci_mlst = create_db("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst")
        download_profiles_st ("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst", loci_mlst)
        print("   ... done \n")
        
        print("Downloading tox database")
        path_tox_sequences, loci_tox = create_db("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        download_profiles_tox ("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        print("   ... done \n")
        
        print ('bash ' + arguments.path + '/data/resistance/update_database_resistance.sh')
        os.system('bash ' + arguments.path + '/data/resistance/update_database_resistance.sh')
        print("   ... done \n\n\n")
    try:
        os.makedirs(arguments.outdir)
        print("Directory '%s' created successfully \n" %arguments.outdir)
    except OSError :
        print("Directory '%s' can not be created \n"  %arguments.outdir)        
        sys.exit(0)