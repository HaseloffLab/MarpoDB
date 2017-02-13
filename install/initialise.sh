#!/bin/bash
# Usage sh initialise.sh [Database_name]
if [ -z "$1" ]
  then
    echo "\n##\n#\n# ERROR - No database name supplied\n#\n# Usage: sh initialise.sh [DATABASE_NAME]\n#\n##\n"
	exit 1;
fi

cd ..

DBname=$1
echo "Loading sequences into database"
python server/initialise.py ${DBname}

echo "Getting sequences from database for annotation"
python server/getsequences.py ${DBname} CDS

mkdir temp
mv CDSs.fa temp/sequences.fa

cd install