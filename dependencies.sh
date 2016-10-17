#!/bin/bash
bash
BASE=$(pwd)
mkdir src
cd src
SRC=$(pwd)

# Transdecoder
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.tar.gz --no-check-certificate
tar xvfz v3.0.0.tar.gz
cd TransDecoder-3.0.0
make
TRANSDECODER=$(pwd)
export PATH="$PATH:$TRANSDECODER/"
echo "Unless errors appeared, TransDecoder successfully installed"
cd $SRC

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

cd $BASE