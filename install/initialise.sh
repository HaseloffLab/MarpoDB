#!/bin/bash
# Usage sh initialise.sh [Database_name]
cd ..

DBname=$1
echo "Loading sequences into database"
python server/initialise.py ${DBname}

echo "Getting sequences from database for annotation"
python server/getsequences.py ${DBname} CDS

mkdir temp
mv CDSs.fa temp/sequences.fa

cd install
