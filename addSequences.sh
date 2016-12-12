#!/bin/bash
# USAGE - sh addSequences.sh [TRANSCRIPTS_FILE] [GENOME_FILE] [DBNAME] [NUM THREADS]

transcriptFile=$1
genomeFile=$2
dbName=$3
numThreads=$4

uniprotFasta=data/Uniprot/uniprot.fasta
pfamHMM=data/Pfam/Pfam-A.hmm

coverageTh=0.99
identityTh=0.99



prefix=${transcriptFile##*/}
gprefix=${genomeFile##*/}

outputDir=output/${prefix}_${gprefix}

mkdir output
mkdir ${outputDir}

# Step 1: CDS prediction using Transdecoder

TransDecoder.LongOrfs -t ${transcriptFile}

makeblastdb -in ${uniprotFasta} -dbtype 'prot'

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

makeblastdb -in data/splign_temp/fa_seq/genome.fa -dbtype nucl

makeblastdb -in data/splign_temp/fa_seq/transcripts.fa -dbtype nucl

compart -qdb data/splign_temp/fa_seq/transcripts.fa -sdb data/splign_temp/fa_seq/genome.fa > data/splign_temp/transcripts.compartments

python scripts/fixCompartmentNames.py data/splign_temp/fa_seq/genome.fa data/splign_temp/fa_seq/transcripts.fa data/splign_temp/transcripts.compartments > data/splign_temp/transcripts.compartments.fixed

splign -ldsdir data/splign_temp/fa_seq -comps data/splign_temp/transcripts.compartments.fixed > ${outputDir}/${prefix}_${gprefix}.splign

# rm -rf data/splign_temp

python scripts/formatSplign.py ${outputDir}/${prefix}_${gprefix}.splign > ${outputDir}/${prefix}_${gprefix}.splign.gff3 ${coverageTh} ${identityTh}

python scripts/filterFastaByMap.py  ${outputDir}/${prefix}.transdecoder.complete.pep ${outputDir}/${prefix}_${gprefix}.splign.gff3 >  ${outputDir}/${prefix}.transdecoder.complete.mapped.pep

python scripts/filterFastaByMap.py  ${outputDir}/${prefix}.transdecoder.complete.trans ${outputDir}/${prefix}_${gprefix}.splign.gff3 >  ${outputDir}/${prefix}.transdecoder.complete.mapped.trans

python scripts/filterPfamByMap.py ${outputDir}/${prefix}.Pfam.domtblout ${outputDir}/${prefix}_${gprefix}.splign.gff3 > ${outputDir}/${prefix}.Pfam.mapped.domtblout

# Step 3: Loading data into DB 

python scripts/cutGenes.py ${outputDir}/${prefix}_${gprefix}.splign.gff3  ${outputDir}/${prefix}.transdecoder.complete.mapped.trans ${genomeFile}  ${outputDir}/${prefix}.transdecoder.complete.mapped.pep ${dbName} > ${outputDir}/log.txt

rm genome.fa.*
