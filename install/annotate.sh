#!/bin/bash

# Usage sh annotate.sh [numThreads] [Coverage_threshold] [Identity_threshold]
#
# We suggest to use a minimum coverage_threshold of 20 and a minimum identity_threshold of 35.
#
cd ..

tempDir=temp
numThreads=$1
cov=$2
id=$3
blastDB=data/Uniprot/uniprot.fasta

if [ -z "$1" ]
  then
    echo "\n##\n#\n# ERROR - No arguments supplied\n#\n# Usage: sh annotate.sh [NUMBER_THREADS] [COVERAGE_THRESHOLD] [IDENTITY_THRESHOLD]\n#\n# We suggest using a minimum coverage_threshold of 20 and a minimum identity_threshold of 35.\n#\n##\n"
	exit 1;
fi

mkdir ${tempDir}

# Blastp search
makeblastdb -in ${blastDB} -dbtype 'prot'

blastp -query ${tempDir}/sequences.fa -db data/Uniprot/uniprot.fasta -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send qcovs pident evalue" -evalue 1e-6 -num_threads ${numThreads} > temp/blastOut.outfmt6

python scripts/filterBlastByCl.py ${tempDir}/blastOut.outfmt6 ${cov} ${id} > ${tempDir}/blastOut.filtered.outfmt6

python scripts/addBlastInfo.py ${tempDir}/blastOut.filtered.outfmt6 data/Uniprot/uniprot.info > ${tempDir}/blastp.info

# Pfam motif search
hmmpress data/Pfam/Pfam-A.hmm

hmmscan --cpu ${numThreads} --domtblout ${tempDir}/Pfam.domtblout --cut_ga data/Pfam/Pfam-A.hmm temp/sequences.fa

# InterproScan
mkdir ${tempDir}/interpro
cat ${tempDir}/sequences.fa | sed 's/*//g' > ${tempDir}/sequences_clean.fa
interproscan.sh -d ${tempDir}/interpro -f gff3 html -goterms -pa -i ${tempDir}/sequences_clean.fa

cd ${tempDir}/interpro
tar xvfz sequences_clean.fa.html.tar.gz
mkdir ../../server/templates/interpro
ls *.html > list
while read line; do ruby ../../scripts/parseInterpro.rb $line > "../../server/templates/interpro/"$line; done < list

mkdir data/filtered
cp temp/Pfam.domtblout data/filtered/Pfam.domtblout
cp temp/blastp.info data/filtered/blastp.info

cd install
