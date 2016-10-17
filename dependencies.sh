#!/bin/bash
bash
mkdir src
cd src

# Transdecoder
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.tar.gz --no-check-certificate
tar xvfz v3.0.0.tar.gz
cd TransDecoder-3.0.0
make
TRANSDECODER=$(pwd)
export PATH="$PATH:$TRANSDECODER/"
echo "Unless errors appeared, TransDecoder successfully installed"
cd ..;

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
cd ..

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
cd ..

# splign and compart
curl ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/CURRENT/ncbi_cxx--12_0_0.tar.gz --user anonymous: -o ncbi_cxx--12_0_0.tar.gz
tar xvfzp ncbi_cxx--12_0_0.tar.gz
cd ncbi_cxx--12_0_0
NCBITOOLS=$(pwd)
./configure --prefix="$NCBITOOLS/"
make
make check
make install
echo "Unless errors appeared, splign and compart should have been installed" 
# InterproScan

curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz -o interproscan-5.20-59.0-64-bit.tar.gz
curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.20-59.0/interproscan-5.20-59.0-64-bit.tar.gz.md5 -o interproscan-5.20-59.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.20-59.0-64-bit.tar.gz.md5
# Must say OK
tar xvzfp interproscan-5.20-59.0-64-bit.tar.gz
# Decide if you want to include Panther models, which requires an additional 15Gb file