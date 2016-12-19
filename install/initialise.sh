#!/bin/bash
# Usage sh initialise.sh [Database_name]
cd ..

DBname=$1
echo "Loading sequences into database"
python install/database/initialise.py ${DBname}

echo "Getting sequences from database for annotation"
python install/database/getsequences.py ${DBname}

cd install