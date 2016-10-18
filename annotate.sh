#!/bin/bash
# Usage sh addSequences.sh [list_file] [dbName] [numThreads]

tempDir=temp
listFile=$1
dbName=$2
numThreads=$3
blastDB=data/Uniprot/uniprot.fasta

mkdir ${tempDir}

python scripts/getSequences.py ${listFile} ${tempDir}/sequences.fa ${dbName} cds
makeblastdb -in ${blastDB} -dbtype 'prot'

blastp -query ${tempDir}/sequences.fa -db data/Uniprot/uniprot.fasta -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send qcovs pident evalue" -evalue 1e-6 -num_threads ${numThreads} > temp/blastOut.outfmt6

python scripts/filterBlastByCl.py ${tempDir}/blastOut.outfmt6 20 35 > ${tempDir}/blastOut.filtered.outfmt6

python scripts/addBlastInfo.py ${tempDir}/blastOut.filtered.outfmt6 data/Uniprot/uniprot.info > ${tempDir}/blastp.info

python scripts/parseBlastp.py ${tempDir}/blastp.info ${dbName}

hmmpress data/Pfam/Pfam-A.hmm

hmmscan --cpu ${numThreads} --domtblout ${tempDir}/Pfam.domtblout --cut_ga data/Pfam/Pfam-A.hmm

python scripts/parsePfam.py ${tempDir}/Pfam.domtblout marpodb