# MarpoDB
An open registry of Marchantia polymorpha genetic parts

# Description
MarpoDB is a gene-centric database developed for genetic engineering and synthetic biology. This is the result of dealing with highly fragmented genomic data (from a non-sequenced organism, Marchantia polymorpha) and compiling it into an accessible resource for sequence exploration and retrieval. The database framework, however, can be used with any type of genetic data and can be set up locally.

In brief, we start off from two "raw" sequence files in FASTA format containing genomic contigs/scaffols and transcripts. Using these files we will perform several analyses, such as ORF prediction, BLAST homology search, HMMR protein motif prediction, mapping transcripts to the "genome" and interrelating resulting data for loading it into a database.

Then, we shall access the database by means of web server accessible by your browser.

# Installation
To install a gene-centric database harboring DNA parts we need to install serveral 3rd party software and databases, libraries for running the web server and application, and finally compile the data before loading into the database. 

Currently, we supply scripts for compiling most libraries for a 64 bit architecture in a linux server, however it may be possible to build the required libraries in a different architecture and perform the corresponding analyses successfully. Also, we are setting up local versions of all the software and libraries, without the need for sudo access since this is the common case for shared servers in academia.

We supply several scripts to split and automate (as much as possible) the setup process.

The setup is divided into the following sections:

## Bioinformatics software and required databases

A wrapper script is provided in dependencies.sh, but if anything fails, try the step-by-step approach.

```bash
# Go to bash first
bash
nohup sh dependencies.sh
```

### Step-by-step approach

Let's set up the directory structure and environmental path variables:
```bash
#!/bin/bash
# Ensuring we are in BASH.
bash

# Setting up environmental variables. Remember to set them again if installing in more than one session and be in the root of the MarpoDB directory. 

# Also if you have already installed everything remember to activate the virtual environment so that the python libraries are referenced correctly: $ source ~/ENV/bin/activate;

BASE=$(pwd)
mkdir src
mkdir data

cd src
SRC=$(pwd)
cd $BASE
cd data
DATA=$(pwd)
cd $BASE

# Adding new custom paths into different files to avoid touching local .profile and duplicating paths. Only run once.
echo "Adding new custom paths into different files to avoid touching local .profile and duplicating paths"
echo "# Adding new custom paths" >> ~/.bashrc
echo ". ~/.paths" >> ~/.bashrc
echo ". ~/.pypaths" >> ~/.bashrc
echo ". ~/.ldpaths" >> ~/.bashrc
```

### Software:

- TransDecoder (part of the Trinity pipeline by Brian Hass -https://github.com/TransDecoder/TransDecoder/releases)
```bash
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
```

- BLAST-2.2.27 (specific version in ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ - https://www.ncbi.nlm.nih.gov/books/NBK279690/)
```bash
## BLAST 2.2.27
cd $SRC
curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ncbi-blast-2.2.27+-x64-linux.tar.gz --user anonymous: -o ncbi-blast-2.2.27+-x64-linux.tar.gz
tar xvfzp ncbi-blast-2.2.27+-x64-linux.tar.gz 
cd ncbi-blast-2.2.27+
BLAST=$(pwd)
export PATH=$BLAST/bin:$PATH
echo "export PATH=$BLAST/bin:$PATH" > ~/.paths
echo "Unless errors appeared, BLAST 2.2.27 successfully installed"

```

- HMMR (http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
```bash
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
```

### Databases:

Now we'll download 3rd party databases and put them in the data/ directory

- Uniprot (TreMBLe and Swissprot - http://www.uniprot.org/downloads)

You might want to include taxonomy filters for obtaining sequences to a specific taxonomy level for you organism and download only those sequences. To do this go to http://www.uniprot.org/uniprot/#orgViewBy. Then, find you taxonomy level and click. Finally choose "Map to - UniProtKB" on the left side panel. Place the file in the data/Uniprot/ directory and rename it to uniprot.fasta .

Otherwise just do:

```bash
## Uniprot - build your own taxonomy database or just download the whole thing (16Gb, 30Gb decompressed for Trembl and 85Mb for Swissprot)
cd $DATA
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -dv uniprot_trembl.fasta.gz
gunzip -dv uniprot_sprot.fasta.gz
cat uniprot_sprot.fasta >> uniprot_trembl.fasta
mv uniprot_trembl.fasta uniprot.fasta
rm uniprot_sprot.fasta
mkdir Uniprot
mv uniprot.fasta Uniprot/
```

- Pfam-A (check most recent release hmm version - http://pfam.xfam.org FTP tab)
```bash
## Pfam
cd $DATA
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
gunzip -dv Pfam-A.hmm.gz
mkdir Pfam
mv Pfam-A.hmm Pfam/
```

- pip (python installer) and setting up virtualenv
```bash
# pip
easy_install --user pip
export PATH=$PATH:~/.local/bin
echo "export PATH=$PATH:~/.local/bin" > ~/.paths
pip install virtualenv --user
virtualenv ~/ENV
```
Activate the virtualenv. This has to be performed in every new session for the system to locate the virtualenv directory ~/ENV

```bash
source ~/ENV/bin/activate
```

- python libraries
```bash
# python libraries to ~/ENV virtualenvironment
cd $BASE
pip install -r requirements.txt

```

In case anything fails, go to the original source and read the documentation. Errors may be due to new releases or compilation errors but should be straightforward to determine.

# Problematic libraries and software

These ones have some issues locating system libraries and might generate errors in older distributions. We provide simple "hacks" to get them up an running but in the case you encounter any errors you might want to google the ways around them...

- PostgreSQL (check most recent release - https://www.postgresql.org/download/ . For this tutorial we'll compile from source)
```bash
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
```
Now lets try to set up the server and create a psql database with appropriate credentials

```bash
## Lest start postgres and make a database and a user
initdb -D ~/var/ -U postgres
pg_ctl -D /disk1/bp358/var/ -l logfile start
# Log into postgres
psql -U postgres
```
Once logged in to psql, create a database with the name desired for loading the data (will be required afterwards) and a user with priviledge to write on it. Some distributions directly refer the logged user as the "role" required. If that is the case then use that user.
```SQL
CREATE DATABASE test;
CREATE ROLE <user>;
ALTER ROLE "<user>" with LOGIN;
\q;
```

# Psycopg python library compilation from source
```bash
# psycopg from source
cd $SRC
wget http://initd.org/psycopg/tarballs/PSYCOPG-2-6/psycopg2-2.6.2.tar.gz
tar xvfzp psycopg2-2.6.2.tar.gz 
cd psycopg2-2.6.2
python setup.py build_ext --pg-config $POSTGRESQL/bin/pg_config  --build-lib $POSTGRESQL/lib build && python setup.py install --user
cd build
PYTHONP=$(pwd)
export PYTHONPATH=${PYTHONPATH}:$PYTHONP
echo "export PYTHONPATH=${PYTHONPATH}:$PYTHONP" > ~/.pypaths
```

- Splign and Compart (NCBI tools - https://www.ncbi.nlm.nih.gov/Web/Newsltr/V14N2/splign.html. We are using a precompiled binary provided from http://sing.citi.uvigo.es/static/BDBM/ncbi.tar.gz)
```bash
cd $SRC
wget http://sing.citi.uvigo.es/static/BDBM/ncbi.tar.gz
tar xfvzp ncbi.tar.gz
mkdir ncbi_bins
mv splign ncbi_bins/
mv compart ncbi_bins/
cd ncbi_bins
BINS=$(pwd)
export PATH=$PATH:$BINS
echo "export PATH=$PATH:$BINS" > ~/.paths
```
If splign complains about not being able to find libpcre.so.0 do:

```bash
#LIBPCRE=$(locate libpcre.so.3 | head -1)		
#cp $LIBPCRE .										
#mv libpcre.so.3 libpcre.so.0					
#export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH
#echo "export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH" ~/.ldpaths
```

Lets check if splign works:
```bash
echo "Let's check everything works. Now will try to run splign..."
splign -help | head -20
echo "Unless errors appeared, Splign and compart successfully installed"
```

- InterproScan (5Gb - Check the requirements in https://github.com/ebi-pf-team/interproscan/wiki)
The lastest version of InterproScan uses a JAVA 1.8, be sure to select the appropriate InterproScan version so that you won't run into any problems running it.

```bash
## InterproScan 5.20-59.0 (Java 8)
cd $SRC
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz -o interproscan-5.20-59.0-64-bit.tar.gz
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5 -o interproscan-5.20-59.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5
## Must say OK, otherwise download again
sleep 5
tar xvzfp interproscan-5.20-59.0-64-bit.tar.gz
cd interproscan-5.20-59.0-64-bit
```

If this fails because you have other version of Java, try:

```bash
## InterproScan 5.16-55.0 (Java 6/7)
cd $SRC
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.16-55.0-64-bit.tar.gz.md5
tar xvfzp interproscan-5.16-55.0-64-bit.tar.gz
cd interproscan-5.16-55.0-64
```

Now, decide if you want Panther models included, since they take quite a lot of space (15Gb)

```bash
## Decide if you want to include Panther models, which requires an additional 15Gb file
cd data
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-10.0.tar.gz.md5
md5sum -c panther-data-10.0.tar.gz.md5
## This must say OK, otherwise download again
tar xfvzp panther-data-10.0.tar.gz
```

Look-up service... This is provided by InterPro as a web service in EBI if you can access or want to have a local look-up table, you'll need 96Gb! otherwise you might want to disable it...

https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload

Now let's check if interproscan works

```bash
echo "Now let's check if InterproScan works"
./interproscan.sh
INTERPRO=$(pwd)
export PATH=$PATH:$INTERPRO
echo "export PATH=$PATH:$INTERPRO" > ~/.paths
```

## Data compilation

cd $BASE
```bash
nohup sh addSequences.sh <TRANSCRIPTS FILE> <GENOME FILE> <DATABASE NAME> <NUMBER OF PROCESSORS> &
```

