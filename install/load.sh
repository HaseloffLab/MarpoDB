#!/bin/bash
# Usage sh load.sh [Database_name]
DBname=$1

if [ -z "$1" ]
  then
    echo "\n##\n#\n# ERROR - No database name supplied\n#\n# Usage: sh load.sh [DATABASE_NAME]\n#\n##\n"
	exit 1;
fi

cd ..

echo "Loading annotations to database"
python server/annotate.py ${DBname}

cd install
