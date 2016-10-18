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

cd $SRC
## TransDecoder
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.tar.gz --no-check-certificate
tar xvfz v3.0.0.tar.gz
cd TransDecoder-3.0.0
make
TRANSDECODER=$(pwd)
echo "export PATH=$PATH:$TRANSDECODER" >> ~/.bashrc
echo "Unless errors appeared, TransDecoder successfully installed"
cd $SRC

## BLAST 2.2.27
curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ncbi-blast-2.2.27+-x64-linux.tar.gz --user anonymous: -o ncbi-blast-2.2.27+-x64-linux.tar.gz
tar xvfzp ncbi-blast-2.2.27+-x64-linux.tar.gz 
cd 
BLAST=$(pwd)
./configure --prefix="$BLAST/"
make
make install
export PATH=$PATH$BLAST/bin
echo "export PATH=$PATH:$BLAST/bin" >> ~/.bashrc
echo "Unless errors appeared, BLAST 2.2.27 successfully installed"
cd $SRC

## HMMR
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar xfvzp hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd hmmer-3.1b2-linux-intel-x86_64
HMMR=$(pwd)
./configure --prefix=$HMMR
make check
make install
export PATH="$PATH:$HMMR/bin"
echo "export PATH=$PATH:$HMMR/bin" >> ~/.bashrc
echo "Unless errors appeared, HMMR successfully installed"
cd $SRC

## Splign and compart
wget http://sing.citi.uvigo.es/static/BDBM/ncbi.tar.gz
tar xfvzp ncbi.tar.gz
mkdir ncbi_bins
mv splign ncbi_bins/
mv compart ncbi_bins/
cd ncbi_bins
BINS=$(pwd)
export PATH=$PATH:$BINS
echo "export PATH=$PATH:$BINS" >> ~/.bashrc
echo "Let's check everything works. Now will try to run splign..."
sleep 5
splign -help | head -20
echo "Unless errors appeared, Splign and compart successfully installed"

## If splign complains about not being able to find libpcre.so.0 do

#LIBPCRE=$(locate libpcre.so.3 | head -1)		
#cp $LIBPCRE .										
#mv libpcre.so.3 libpcre.so.0					
#export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH
#echo "export LD_LIBRARY_PATH=$BINS:$LD_LIBRARY_PATH" >> ~/.bashrc

## InterproScan 5.20-59.0 (Java 8)
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz -o interproscan-5.20-59.0-64-bit.tar.gz
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5 -o interproscan-5.20-59.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5
## Must say OK, otherwise download again
sleep 5
tar xvzfp interproscan-5.20-59.0-64-bit.tar.gz
cd interproscan-5.20-59.0-64-bit

## InterproScan 5.16-55.0 (Java 6/7)
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.16-55.0/interproscan-5.16-55.0-64-bit.tar.gz.md5
#md5sum -c interproscan-5.16-55.0-64-bit.tar.gz.md5
#tar xvfzp interproscan-5.16-55.0-64-bit.tar.gz
#cd interproscan-5.16-55.0-64-bit

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
echo "export PATH=$PATH:$INTERPRO" >> ~/.bashrc

## Databases

cd $DATA

## Uniprot - build your own taxonomy database or just download the whole thing (16Gb, 30Gb decompressed for Trembl and 85Mb for Swissprot)
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -dv uniprot_trembl.fasta.gz
gunzip -dv uniprot_sprot.fasta.gz
cat uniprot_sprot.fasta >> uniprot_trembl.fasta
mv uniprot_trembl.fasta uniprot.fasta
rm uniprot_sprot.fasta
mkdir Uniprot
mv uniprot.fasta Uniprot/
cd $DATA

## Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm.gz
gunzip -dv Pfam-A.hmm.gz
mkdir Pfam
mv Pfam-A.hmm Pfam/