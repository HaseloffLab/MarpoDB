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

# Bioinformatics software and libraries

We need to install the following software:

- TransDecoder (part of the Trinity pipeline by Brian Hass -https://github.com/TransDecoder/TransDecoder/releases)

- BLAST-2.2.27 (specific version in ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.27/ - https://www.ncbi.nlm.nih.gov/books/NBK279690/)

- HMMR (http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)

- Splign (NCBI tools - https://www.ncbi.nlm.nih.gov/Web/Newsltr/V14N2/splign.html)

- InterproScan (5Gb - Check the requirements in https://github.com/ebi-pf-team/interproscan/wiki)

We provide a script for installing into a src/ directory. In case anything fails, just read the dependencies.sh file and check which step failed.

```bash
sh dependencies.sh
```

# Data compilation




# Server

python2
PostgreSQL
Flask
Jinja2
Biopython


