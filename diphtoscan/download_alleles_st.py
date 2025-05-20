#!/usr/bin/env python3
# Example script to download alleles from a sequence definition database
# Written by Keith Jolley
# Copyright (c) 2017, University of Oxford
# E-mail: keith.jolley@zoo.ox.ac.uk
#
# This file is part of Bacterial Isolate Genome Sequence Database (BIGSdb).
#
# BIGSdb is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BIGSdb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import os
import requests
import pandas as pd
import io 

BASE_URI = 'https://bigsdb.pasteur.fr/api'

def download_alleles(database:str, scheme_id:str, folder:str) -> list:
    """
    Parameters
    ----------
    database  : Database configuration name
    scheme_id : Only return loci belonging to scheme. If this option is 
                not used then all loci from the database will be downloaded
    folder    : Output directory

    Returns
    -------
    Name of downloaded loci 
    """
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
    dir = folder or './'
    url = BASE_URI + '/db/' + database
    r = requests.get(url)
    if r.status_code == 404:
        print('Database ' + database + ' does not exist.')
        os._exit(1)
    loci = []
    if scheme_id:
        url = BASE_URI +  '/db/' + database + '/schemes/' + str(scheme_id);
        r = requests.get(url);
        if r.status_code == 404:
            print('Scheme ' + str(scheme_id) + ' does not exist.');
            os._exit(1)
        loci = r.json()['loci']
    else:
        url = BASE_URI + '/db/' + database + '/loci?return_all=1'
        r = requests.get(url);
        loci = r.json()['loci'];
    name_loci = []
    for locus_path in loci:
        r = requests.get(locus_path)
        locus = r.json()['id']
        if r.json()['alleles_fasta']:
            r = requests.get(r.json()['alleles_fasta'])
            name_loci.append(locus)
            fasta_file = open(dir + '/' + locus + '.fas', 'w')
            fasta_file.write(r.text)
            fasta_file.close()
    return name_loci


def create_db (database:str, scheme_id:str, folder:str):
    loci_mlst = download_alleles(database, scheme_id, folder+"/sequences")
    path_loci_mlst = [folder+"/sequences/"+ locus +'.fas' for locus in loci_mlst]
    path_database = folder +"/"+ database +"_scheme_"+ scheme_id+ ".fas"
    os.system("cat "+" ".join(path_loci_mlst)+" > "+ path_database + " 2>/dev/null")
    return path_database, loci_mlst


def download_profiles_st (database:str, scheme_id:str, folder:str, loci_mlst:list):
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
    dir = folder or './'
    url = BASE_URI + '/db/' + database
    r = requests.get(url)
    if r.status_code == 404:
        print('Database ' + database + ' does not exist.')
        os._exit(1)
    if scheme_id:
        url = BASE_URI +  '/db/' + database + '/schemes/' + str(scheme_id);
        r = requests.get(url+"/profiles_csv");
        table_profiles_st = pd.read_csv(io.StringIO(r.text), sep="\t", index_col=0, dtype=str)
        table_profiles_st[loci_mlst].to_csv(dir + '/st_profiles.txt', sep='\t')  
    return dir + '/st_profiles.txt'

def download_profiles_tox (database:str, scheme_id:str, folder:str):
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
    dir = folder or './'
    url = BASE_URI + '/db/' + database
    r = requests.get(url)
    if r.status_code == 404:
        print('Database ' + database + ' does not exist.')
        os._exit(1)
    if scheme_id:
        url = BASE_URI +  '/db/' + database + '/schemes/' + str(scheme_id);
        r = requests.get(url+"/profiles_csv");
        table_profiles_st = pd.read_csv(io.StringIO(r.text), sep="\t", index_col=0, dtype=str)
        table_profiles_st['tox'].to_csv(dir + '/tox_profiles.txt', sep='\t')  
    return dir + '/tox_profiles.txt'