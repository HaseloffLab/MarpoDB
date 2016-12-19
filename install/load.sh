#!/bin/bash
# Usage sh load.sh [Database_name]
DBname=$1

cd ..

echo "Loading annotations to database"
python install/database/annotate.py ${DBname}

cd install
