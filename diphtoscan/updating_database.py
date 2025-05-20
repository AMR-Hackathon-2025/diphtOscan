import datetime
import re
import subprocess

from pathlib import Path

import requests

import pandas as pd

from .download_alleles_st import create_db, download_profiles_st, download_profiles_tox

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


def download_amrfinder_database(url: str, path: str) -> None:
    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    href_re = re.compile(r'<a href="[^"]+">([^<]+)</a>')
    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download file from {url}: {response.status_code}")
    start_parsing = False
    for line in response.iter_lines():
        if line.startswith(b'<pre>Name'):
            start_parsing = True
        elif line.startswith(b'<hr></pre>'):
            break
        if start_parsing and line.startswith(b'<a href='):
            match = href_re.search(line.decode('utf-8'))
            if match:
                filename = match.group(1)
                if not filename.endswith('/'):
                    file_url = url + filename
                    response = requests.get(file_url, stream=True)
                    if response.status_code == 200:
                        with open(output_dir / filename, 'wb') as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                f.write(chunk)
                    else:
                        raise RuntimeError(f"Failed to download file from {file_url}: {response.status_code}")


def update_amrfinderplus_db_file(input_path: str, output_path: str, skip_first_line: bool = False) -> None:
    """
    Update the AMRFinderPlus database file by appending lines from the input file to the output file.
    """
    with open(input_path) as input_file:
        with open(output_path, 'a') as output_file:
            if skip_first_line:
                next(input_file)  # Skip the first line
            for line in input_file:
                output_file.write(line)



def remove_mlst_database(mlst_db_path: str):
    """
    Remove the MLST database files from the specified path.
    """
    db_path = Path(mlst_db_path).parent
    for path in db_path.glob('*'):
        if not path.is_dir():
            path.unlink()
    for path in (db_path / 'sequences').glob('*'):
        path.unlink()
    (db_path / 'sequences').rmdir()


def update_database(arguments, mlst_database: tuple[str, str], tox_database:tuple):
    if arguments.update :
        date = datetime.datetime.today().strftime('%Y-%m-%d') 

        # remove old database
        remove_mlst_database(mlst_database[1])
        print("Downloading MLST database")
        path_mlst_sequences, loci_mlst = create_db("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst")
        download_profiles_st ("pubmlst_diphtheria_seqdef", "3", arguments.path +"/data/mlst", loci_mlst)
        print("   ... done \n")

        remove_mlst_database(tox_database[1])
        print("Downloading tox database")
        path_tox_sequences, loci_tox = create_db("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        download_profiles_tox ("pubmlst_diphtheria_seqdef", "4", arguments.path +"/data/tox")
        print("   ... done \n")
        
        # wget --quiet --recursive --no-parent --no-host-directories --cut-dirs=6 -e robots=off https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ -P $PATH_DB/$DATE 

        # find AMRFinderPlus version
        amrfinderplus_version = subprocess.run(['amrfinder', '--version'], capture_output=True, text=True).stdout.split('.')[0]
        if amrfinderplus_version == '3':
            # URL of latest AMRFinderPlus 3 compatible database
            url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-07-22.1/'
            tsv_suffix = 'tab'
        elif amrfinderplus_version == '4':
            url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/'
            tsv_suffix = 'tsv'
        else:
            raise RuntimeError(f"Unsupported AMRFinderPlus version: {amrfinderplus_version}")
        
        print(f"Downloading AMRFinderPlus database and saving to {amr_database_path}")
        source_amr_database_path = arguments.path + '/data/resistance/Corynebacterium_diphtheriae'
        amr_database_path = arguments.path + '/data/resistance/' + date
        download_amrfinder_database(url, amr_database_path)
        open(amr_database_path + '/version.txt', 'w').write(date + '.1')

        print("Merging custom database with AMRFinderPlus database")
        # conctatenate proteins to AMRFinderPlus database
        fam_file_path = amr_database_path + '/fam.' + tsv_suffix
        update_amrfinderplus_db_file(source_amr_database_path + '/AMRProt_Cd', amr_database_path + '/AMRProt.fa')
        update_amrfinderplus_db_file(source_amr_database_path + '/fam_Cd.tab',
                                     fam_file_path, skip_first_line=True)
        complete_missing_classification(fam_file_path)

        print("Building BLAST database for AMRFinderPlus")
        # makeblastdb -in $PATH_DB/$DATE/AMRProt -dbtype prot  -logfile /dev/null
        subprocess.run(['makeblastdb', '-in', amr_database_path + '/AMRProt.fa',
                        '-dbtype', 'prot', '-out', amr_database_path + '/AMRProt',
                        '-logfile', '/dev/null'])
        print("   ... done \n\n\n")