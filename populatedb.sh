# A script for populating psql database. Create a database by running createdb marpodb, then create tables by running psql -d marpodb -f ../psql/tables.sql
psql -d marpodb -f psql/tables.sql
python cutGenes.py data/Cam1_fpkm0.1_iso10.trans.mRNA_to_Cam1_Merac.splign.gff3 ../data/Cam1_Br_fpkm0.1_isopct10_transcripts.fasta ../data/final.scaffolds.fa ../data/Cam1_Br_filtered_fpkm0.1_isopct10.fasta.transdecoder.pep
python parseFasta.py data/Cam1_Br_filtered_fpkm0.1_isopct10.fasta.transdecoder.pep
python parsePfam.py data/pfam_ga.dombtlout
python parseBlastp.py data/blastp.info
python parseBlastm.py data/blastm.info
