#!/bin/bash

PATH_DB=$(dirname "$0")
DATE=$(date "+%Y-%m-%d")

echo "Indexing" ;

hmmpress -f $PATH_DB/$DATE/AMR.LIB > /dev/null 2> /dev/null

makeblastdb -in $PATH_DB/$DATE/AMRProt -dbtype prot  -logfile /dev/null
makeblastdb -in $PATH_DB/$DATE/AMR_CDS -dbtype nucl  -logfile /dev/null

taxgroups=$(awk '{if ($3>0 && $1!="#taxgroup") print $1}' $PATH_DB/$DATE/taxgroup.tab)
for taxgroup in $taxgroups  
do makeblastdb -in $PATH_DB/$DATE/AMR_DNA-$taxgroup -dbtype nucl  -logfile /dev/null
done

echo -e "Corynebacterium_diphtheriae\tCorynebacterium_diphtheriae\t0" >> $PATH_DB/$DATE/taxgroup.tab

PATH_DB="$PATH_DB/$DATE" 
VERSION="$DATE"
echo "Database directory: '$PATH_DB'"
echo "Database version: $DATE.1"