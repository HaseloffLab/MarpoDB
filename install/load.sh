#!/bin/bash
# Usage sh load.sh [Database_name]
DBname=$1

cd ..

echo "Loading annotations to database"
python server/annotate.py ${DBname}

cd install
