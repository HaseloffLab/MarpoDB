# MarpoDB
An open registry of Marchantia polymorpha genetic parts

# Description
MarpoDB is a gene-centric database developed for genetic engineering and synthetic biology. This is the result of dealing with highly fragmented genomic data (from a non-sequenced organism, Marchantia polymorpha) and compiling it into an accessible resource for sequence exploration and retrieval. The database framework, however, can be used with any type of genetic data and can be set up locally.

In brief, we start off from two "raw" sequence files in FASTA format containing genomic contigs/scaffols and transcripts. Using these files we will perform several analyses, such as ORF prediction, BLAST homology search, HMMR protein motif prediction, mapping transcripts to the "genome" and interrelating resulting data for loading it into a database.

Then, we shall access the database by means of web server accessible by your browser.

# Installation
To install a gene-centric database harboring DNA parts we need to install serveral 3rd party libraries, software and compile the data before loading into the database. 

Currently, we supply scripts for compiling the libraries in a 64 bit architecture run in a linux server, however it may be possible to build the required libraries in a different architecture and perform the corresponding analyses successfully.

For this purpose we supply several scripts to split and automate (as much as possible) the setup process.

The setup is divided into the following sections:

## Bioinformatics software and required databases

A wrapper script is provided in dependencies.sh, but if anything fails, try the step-by-step approach.

```bash
sh dependencies.sh
```

* Step-by-step

Let's set up the directory structure and environmental path variables:
```bash
#!/bin/bash
bash
BASE=$(pwd)
mkdir src
cd src
SRC=$(pwd)
cd $BASE
mkdir data
cd data
DATA=$(pwd)
cd $BASE
```

### Software:

- TransDecoder (part of the Trinity pipeline by Brian Hass -https://github.com/TransDecoder/TransDecoder/releases)
```bash
# Transdecoder
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.tar.gz --no-check-certificate
tar xvfz v3.0.0.tar.gz
cd TransDecoder-3.0.0
make
TRANSDECODER=$(pwd)
export PATH="$PATH:$TRANSDECODER/"
echo "Unless errors appeared, TransDecoder successfully installed"
cd $SRC
```

- BLAST-2.2.27 (specific version in ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ - https://www.ncbi.nlm.nih.gov/books/NBK279690/)
```bash
# BLAST 2.2.27
curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ncbi-blast-2.2.27+-x64-linux.tar.gz --user anonymous: -o ncbi-blast-2.2.27+-x64-linux.tar.gz
tar xvfzp ncbi-blast-2.2.27+-x64-linux.tar.gz 
cd 
BLAST=$(pwd)
./configure --prefix="$BLAST/"
make
make install
export PATH="$PATH:$BLAST/bin"
echo "Unless errors appeared, BLAST 2.2.27 successfully installed"
cd $SRC
```

- HMMR (http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
```bash
# HMMR
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar xfvz 
cd hmmer-3.1b2-linux-intel-x86_64
HMMR=$(pwd)
./configure --prefix=$HMMR
make check
make install
export PATH="$PATH:$HMMR/bin"
echo "Unless errors appeared, HMMR successfully installed"
cd $SRC
```

- Splign and Compart (NCBI tools - https://www.ncbi.nlm.nih.gov/Web/Newsltr/V14N2/splign.html. We are using a precompiled binary provided from http://sing.citi.uvigo.es/static/BDBM/ncbi.tar.gz)
```bash
# Splign and compart
wget http://sing.citi.uvigo.es/static/BDBM/ncbi.tar.gz
tar xfvzp ncbi.tar.gz
mkdir ncbi_bins
mv splign ncbi_bins/
mv compart ncbi_bins/
cd ncbi_bins
BINS=$(pwd)
export PATH="$PATH:$BINS"
echo "Let's check everything works. Now will try to run splign..."
sleep 5
splign -help | head -20
echo "Unless errors appeared, Splign and compart successfully installed"

## If splign complains about not being able to find libpcre.so.0 do
#echo "LIBPCRE=$(locate libpcre.so.3 | head -1)		
#cp $LIBPCRE .										
#mv libpcre.so.3 libpcre.so.0					
#export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH
```

- InterproScan (5Gb - Check the requirements in https://github.com/ebi-pf-team/interproscan/wiki)
```bash
# InterproScan
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz -o interproscan-5.20-59.0-64-bit.tar.gz
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5 -o interproscan-5.20-59.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5
# Must say OK, otherwise download again
sleep 5
tar xvzfp interproscan-5.20-59.0-64-bit.tar.gz
# Decide if you want to include Panther models, which requires an additional 15Gb file
cd interproscan-5.20-59.0-64-bit
./interproscan
```

### Databases:

- Uniprot (TreMBLe and Swissprot)

- Pfam-A (check most recent release hmm version)

- nr

```bash
cd $DATA
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
```

In case anything fails, go to the original source and read the documentation. Errors may be due to new releases or compilation errors but should be straightforward to determine.

```bash
sh dependencies.sh
```

# Data compilation




# Server installation

You should have pip and python2 installed for performing automated installation.
```bash
sudo apt-get pip python2 postgresql-9.4
# or use the appropriate package manager
pip install -r requirements.txt

```

