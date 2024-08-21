import datetime
import os
import sys

import pandas as pd

from module.download_alleles_st import create_db, download_profiles_st, download_profiles_tox

node_class = {'pld':'OTHER_TOXINS',
'spaA' : 'SpaA-type_pili_diphtheriae',
'spaB' : 'SpaA-type_pili_diphtheriae',
'spaC' : 'SpaA-type_pili_diphtheriae',
'srtA' : 'SpaA-type_pili_diphtheriae',
'spaD' : 'SpaD-type_pili_diphtheriae',
'spaE' : 'SpaD-type_pili_diphtheriae',
'spaF' : 'SpaD-type_pili_diphtheriae',
'srtB' : 'SpaD-type_pili_diphtheriae',
'srtC' : 'SpaD-type_pili_diphtheriae',
'spaG' : 'SpaH-type_pili_diphtheriae',
'spaH' : 'SpaH-type_pili_diphtheriae',
'spaI' : 'SpaH-type_pili_diphtheriae',
'srtD' : 'SpaH-type_pili_diphtheriae',
'srtE' : 'SpaH-type_pili_diphtheriae',
'tox' : 'TOXIN',
'cbpA' : 'VIRULENCE/ADHESIN',
'nanH' : 'VIRULENCE/ADHESIN',
}

def complete_missing_classification(path:str):
    df = pd.read_csv(path, sep="\t", escapechar="\\", engine="python")
    missing_class = df.loc[df['parent_node_id']=='VIRULENCE_Cdiphth']
    for index in missing_class.index:
        for field in ['class','subclass']:
            if pd.isna(df.iloc[index, df.columns.get_loc(field)]) :
                df.iloc[index, df.columns.get_loc(field)] = node_class[df.iloc[index]['#node_id']]
    df.to_csv(path, sep="\t", escapechar="\\", index=False)
    return


def update_database(arguments, mlst_database:tuple, tox_database:tuple):
    if arguments.update :
        date = datetime.datetime.today().strftime('%Y-%m-%d') 

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
        complete_missing_classification(arguments.path + '/data/resistance/' + date + '/fam.tab')
        os.system('bash ' + arguments.path + '/data/resistance/making_blastdb.sh')
        print("   ... done \n\n\n")