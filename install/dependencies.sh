#!/bin/bash
#
# USAGE - sh dependencies.sh [DATABASE_NAME]
#

if [ -z "$1" ]
  then
    echo "\n##\n#\n# ERROR - No arguments supplied\n#\n# USAGE - sh dependencies.sh [DATABASE_NAME]\n#\n##\n"
	exit 1;
fi

# Ensuring we are in BASH.
bash
cd ..
$DATABASENAME=$1

## SOFTWARE

## Python libraries
cd $BASE
pip install -r requirements.txt

## TransDecoder
cd $SRC
wget https://github.com/TransDecoder/TransDecoder/archive/v2.1.0.tar.gz --no-check-certificate
tar xvfz v2.1.0.tar.gz
cd TransDecoder-2.1.0
make
TRANSDECODER=$(pwd)
export PATH=$PATH:$TRANSDECODER
echo "export PATH=$PATH:$TRANSDECODER" > ~/.paths
echo "Unless errors appeared, TransDecoder successfully installed"

## BLAST 2.2.27
cd $SRC
curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ncbi-blast-2.2.27+-x64-linux.tar.gz --user anonymous: -o ncbi-blast-2.2.27+-x64-linux.tar.gz
tar xvfzp ncbi-blast-2.2.27+-x64-linux.tar.gz 
cd ncbi-blast-2.2.27+
BLAST=$(pwd)
export PATH=$BLAST/bin:$PATH
echo "export PATH=$BLAST/bin:$PATH" > ~/.paths
echo "Unless errors appeared, BLAST 2.2.27 successfully installed"

## HMMR
cd $SRC
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar xfvzp hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64
HMMR=$(pwd)
./configure --prefix=$HMMR
make check
make install
export PATH=$PATH:$HMMR/bin
echo "export PATH=$PATH:$HMMR/bin" > ~/.paths
echo "Unless errors appeared, HMMR successfully installed"

## DATABASES

## Uniprot - build your own taxonomy database or just download the whole thing (16Gb, 30Gb decompressed for Trembl and 85Mb for Swissprot)
## Here we will use a taxonomy filter for "Viridiplantae [33090]" for both the sequence database and the annotations file for that database
cd $DATA
mkdir Uniprot
cd Uniprot
wget "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=taxonomy:%22Viridiplantae%20[33090]%22&fil=&format=fasta&force=yes" -O uniprot.fasta.gz 
gunzip uniprot.fasta.gz

# Remember to use the appropiate annotation database, filtered by taxonomy.
wget "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=taxonomy:%22Viridiplantae%20[33090]%22%20AND%20existence:%22evidence%20at%20protein%20level%22%20OR%20existence:%22evidence%20at%20transcript%20level%22&fil=&format=tab&force=yes&columns=id,protein%20names,genes(PREFERRED),genes(ALTERNATIVE),genes(OLN),organism" -O uniprot.info

## Pfam
cd $DATA
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
gunzip -dv Pfam-A.hmm.gz
mkdir Pfam
mv Pfam-A.hmm Pfam/

## Problematic ones...

# PostgreSQL 
cd $SRC
wget https://ftp.postgresql.org/pub/source/v9.6.0/postgresql-9.6.0.tar.gz --no-check-certificate
wget https://ftp.postgresql.org/pub/source/v9.6.0/postgresql-9.6.0.tar.gz.md5 --no-check-certificate
md5sum -c postgresql-9.6.0.tar.gz.md5 
## This must say OK, otherwise download again
tar zvfxp postgresql-9.6.0.tar.gz
cd postgresql-9.6.0
POSTGRESQL=$(pwd)
./configure --prefix=$POSTGRESQL
make
make install
export PATH=$PATH:$POSTGRESQL/bin
echo "export PATH=$PATH:$POSTGRESQL/bin" > ~/.paths 
export LD_LIBRARY_PATH=$POSTGRESQL/lib:$LD_LIBRARY_PATH
echo "export LD_LIBRARY_PATH=$POSTGRESQL/lib:$LD_LIBRARY_PATH" > ~/.ldpaths

## Lest start postgres and make a database
initdb -D ~/var/ -U postgres
pg_ctl -D ~/var/ -l logfile start
# Log into postgres and create DB
psql -U postgres -v v1=$USER  -v v2=$DATABASENAME -f psql/credentials.sql

# psycopg from source
cd $SRC
wget http://initd.org/psycopg/tarballs/PSYCOPG-2-6/psycopg2-2.6.2.tar.gz
tar xvfzp psycopg2-2.6.2.tar.gz 
cd psycopg2-2.6.2
python setup.py build_ext --pg-config $POSTGRESQL/bin/pg_config  --build-lib $POSTGRESQL/lib build && python setup.py install
cd build
PYTHONP=$(pwd)
export PYTHONPATH=${PYTHONPATH}:$PYTHONP
echo "export PYTHONPATH=${PYTHONPATH}:$PYTHONP" > ~/.pypaths

# PartsDB
cd $SRC
git clone https://github.com/HaseloffLab/PartsDB
cd PartsDB
pip install .

### COMPLICATED ONES

## Splign and compart
#curl ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/splign/linux-i64/splign.tar.gz --user anonymous: -o splign.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/splign/linux-i64/splign.tar.gz --no-check-certificate
mkdir ncbi_bins
mv splign.tar.gz ncbi_bins
cd ncbi_bins
tar xfvzp splign.tar.gz
BINS=$(pwd)
export PATH=$PATH:$BINS
echo "export PATH=$PATH:$BINS" > ~/.paths 

## If splign complains about not being able to find libpcre.so.0 do

#LIBPCRE=$(locate libpcre.so | head -1)		
#cp $LIBPCRE .
#LIBPCRE=$(ls libpcre.*)
#mv $LIBPCRE libpcre.so.0	
#export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH
#echo "export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH" > ~/.ldpaths

echo "Let's check everything works. Now will try to run splign..."
sleep 5
splign -help | head -20
echo "Unless errors appeared, Splign and compart successfully installed"

## InterproScan 5.20-59.0 (Java 8)
cd $SRC
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz -o interproscan-5.20-59.0-64-bit.tar.gz
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5 -o interproscan-5.20-59.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5
## Must say OK, otherwise download again
sleep 5
tar xvzfp interproscan-5.20-59.0-64-bit.tar.gz
cd interproscan-5.20-59.0-64-bit

## InterproScan 5.16-55.0 (Java 6/7)
#cd $SRC
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz.md5
#md5sum -c interproscan-5.16-55.0-64-bit.tar.gz.md5
#tar xvfzp interproscan-5.16-55.0-64-bit.tar.gz
#cd interproscan-5.16-55.0

## Decide if you want to include Panther models, which requires an additional 15Gb file

#cd data
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz.md5
#md5sum -c panther-data-10.0.tar.gz.md5
## This must say OK, otherwise download again
#tar xfvzp panther-data-10.0.tar.gz

## Look-up service... This is provided by InterPro as a web service in EBI if you can access or want to have a local look-up table, you'll need 96Gb! otherwise you might want to disable it...
## https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload

echo "Now let's check if InterproScan works"
./interproscan.sh
INTERPRO=$(pwd)
export PATH=$PATH:$INTERPRO
echo "export PATH=$PATH:$INTERPRO" > ~/.paths
cd $BASE
cd install