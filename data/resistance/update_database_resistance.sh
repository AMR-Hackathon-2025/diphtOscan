#!/bin/bash

PATH_DB=$(dirname "$0")
DATE=$(date "+%Y-%m-%d")
## Update Antimicrobial_resistance/AMRFinderPlus

wget --quiet --recursive --no-parent --no-host-directories --cut-dirs=6 -e robots=off https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ -P $PATH_DB/$DATE 
version_DB=$(cat $PATH_DB/$DATE/version.txt)
echo "Looking up databases at https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/"
echo "Downloading AMRFinder database version $version_DB"

## Update with database Corynebacterium_diphtheriae
echo "Updating AMRFinder database by adding the Corynebacterium_diphtheriae database"
echo $DATE.1 > version.txt
mv version.txt $PATH_DB/$DATE/


cat $PATH_DB/Corynebacterium_diphtheriae/AMRProt_Cd >> $PATH_DB/$DATE/AMRProt 
sed '1d' $PATH_DB/Corynebacterium_diphtheriae/AMRProt-mutation_Cd.tab >> $PATH_DB/$DATE/AMRProt-mutation.tab
#sed '1d' $PATH_DB/Corynebacterium_diphtheriae/AMRProt-susceptible_Cd.tab >> $PATH_DB/$DATE/AMRProt-susceptible.tab
sed '1d' $PATH_DB/Corynebacterium_diphtheriae/fam_Cd.tab >> $PATH_DB/$DATE/fam.tab