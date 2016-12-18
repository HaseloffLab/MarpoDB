#!/bin/bash
# Usage sh addSequences.sh [numThreads] [Coverage_threshold] [Identity_threshold]
#
# We suggest to use a minimum coverage_threshold of 20 and a minimum identity_threshold of 35.
#
cd ..

tempDir=temp
numThreads=$1
cov=$2
id=$3
blastDB=data/Uniprot/uniprot.fasta

mkdir ${tempDir}

blastp -query ${tempDir}/sequences.fa -db data/Uniprot/uniprot.fasta -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send qcovs pident evalue" -evalue 1e-6 -num_threads ${numThreads} > temp/blastOut.outfmt6

python scripts/filterBlastByCl.py ${tempDir}/blastOut.outfmt6 ${cov} ${id} > ${tempDir}/blastOut.filtered.outfmt6

python scripts/addBlastInfo.py ${tempDir}/blastOut.filtered.outfmt6 data/Uniprot/uniprot.info > ${tempDir}/blastp.info

## Might want to parametrise taxonomy!!!
wget "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=taxonomy:%22Viridiplantae%20[33090]%22%20AND%20existence:%22evidence%20at%20protein%20level%22%20OR%20existence:%22evidence%20at%20transcript%20level%22&fil=&format=tab&force=yes&columns=id,protein%20names,genes(PREFERRED),genes(ALTERNATIVE),genes(OLN),organism" -O uniprot.info

hmmpress data/Pfam/Pfam-A.hmm

hmmscan --cpu ${numThreads} --domtblout ${tempDir}/Pfam.domtblout --cut_ga data/Pfam/Pfam-A.hmm temp/sequences.fa

mkdir data/filtered
cp temp/Pfam.domtblout data/filtered/Pfam.domtblout
cp temp/blastOut.filtered.outfmt6 data/filtered/blastp_filtered.outfmt6

cd install