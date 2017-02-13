#!/bin/bash
# Usage sh map.sh [numThreads]
cd ..

if [ -z "$1" ]
  then
    echo "\n##\n#\n# ERROR - No argument supplied\n#\n# Usage: sh map.sh [NUM_THREADS]\n#\n##\n"
	exit 1;
fi

transcriptFile=data/transcripts.fa
genomeFile=data/genome.fa
numThreads=$1

uniprotFasta=data/Uniprot/uniprot.fasta
pfamHMM=data/Pfam/Pfam-A.hmm

coverageTh=0.95
identityTh=0.95

prefix=${transcriptFile##*/}
gprefix=${genomeFile##*/}

outputDir=output/${prefix}_${gprefix}

mkdir output
mkdir ${outputDir}

# Step 1: CDS prediction using Transdecoder

#TransDecoder.LongOrfs -t ${transcriptFile}

#makeblastdb -in ${uniprotFasta} -dbtype 'prot'

blastp -query ${prefix}.transdecoder_dir/longest_orfs.pep -db ${uniprotFasta} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads ${numThreads} > ${outputDir}/${prefix}.TransdecoderBLAST.outfmt6

hmmpress ${pfamHMM}

hmmscan --cpu ${numThreads} --domtblout ${outputDir}/${prefix}.Pfam.domtblout --cut_ga ${pfamHMM} ${prefix}.transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t ${transcriptFile} --retain_pfam_hits ${outputDir}/${prefix}.Pfam.domtblout  --retain_blastp_hits ${outputDir}/${prefix}.TransdecoderBLAST.outfmt6

mv ${prefix}.transdecoder_dir ${outputDir}
mv ${prefix}.transdecoder.* ${outputDir}

python scripts/filterComplete.py ${outputDir}/${prefix} ${transcriptFile}

# Step 2: Mapping

rm -rf data/splign_temp
mkdir data/splign_temp

cp ${genomeFile} data/splign_temp/genome.fa

cp  ${outputDir}/${prefix}.transdecoder.complete.trans data/splign_temp/transcripts.fa

mkdir data/splign_temp/fa_seq

python scripts/prepareFASTAforSplign.py data/splign_temp/genome.fa > data/splign_temp/fa_seq/genome.fa

python scripts/prepareFASTAforSplign.py data/splign_temp/transcripts.fa > data/splign_temp/fa_seq/transcripts.fa

splign -mklds data/splign_temp/fa_seq

formatdb -i data/splign_temp/fa_seq/genome.fa -p F -o T
formatdb -i data/splign_temp/fa_seq/transcripts.fa -p F -o T

compart -qdb data/splign_temp/fa_seq/transcripts.fa -sdb data/splign_temp/fa_seq/genome.fa > data/splign_temp/transcripts.compartments

splign -ldsdir data/splign_temp/fa_seq -comps data/splign_temp/transcripts.compartments > ${outputDir}/${prefix}_${gprefix}.splign

python scripts/formatSplign.py ${outputDir}/${prefix}_${gprefix}.splign ${coverageTh} ${identityTh} > ${outputDir}/${prefix}_${gprefix}.splign.gff3 

# Here we need to change a parameter in filterFastaByMap for it to work.

python scripts/filterFastaByMap.py  ${outputDir}/${prefix}.transdecoder.complete.pep ${outputDir}/${prefix}_${gprefix}.splign.gff3 >  ${outputDir}/${prefix}.transdecoder.complete.mapped.pep

python scripts/filterFastaByMap.py  ${outputDir}/${prefix}.transdecoder.complete.trans ${outputDir}/${prefix}_${gprefix}.splign.gff3 >  ${outputDir}/${prefix}.transdecoder.complete.mapped.trans

mkdir data/mapped

cp ${genomeFile} data/mapped/genome.fa
cp ${outputDir}/${prefix}.transdecoder.complete.mapped.pep data/mapped/pep.fa
cp ${outputDir}/${prefix}.transdecoder.complete.mapped.trans data/mapped/trans.fa
cp ${outputDir}/${prefix}_${gprefix}.splign.gff3 data/mapped/map.gff3

cd install
