import sys
import psycopg2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

dnaFile = sys.argv[1]
pepFile = sys.argv[2]

conn = psycopg2.connect("dbname=marpodb")
cur = conn.cursor()

cur.execute("SELECT name, seq FROM gene")
genes = cur.fetchall()

dnaRecords = []

for gene in genes:
	if len(gene[1]) == 0:
		print gene[0], 'is empty'
	else:
		dnaRecords.append(SeqRecord(Seq(gene[1],IUPAC.unambiguous_dna),id=gene[0], description=""))



cur.execute("SELECT name, seq FROM cds")
peps = cur.fetchall()

pepRecords = []

for pep in peps:
	if len(pep[1]) == 0:
                print pep[0], 'is empty'
        else:
		pepRecords.append(SeqRecord(Seq(pep[1],IUPAC.protein),id=pep[0], description=""))


SeqIO.write(dnaRecords, dnaFile, "fasta")
SeqIO.write(pepRecords, pepFile, "fasta")
