import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

import psycopg2
import os

import subprocess

listFile = open(sys.argv[1])
outFile  = open(sys.argv[2], 'w')
dbName = sys.argv[3]
table = sys.argv[4]

conn = psycopg2.connect("dbname={0}".format(dbName))
cur  = conn.cursor()

records = []

for line in listFile:
	cdsID = line.split()[0]
	queryString = "SELECT id, seq FROM {0} WHERE id = '{1}'".format(table, cdsID)
	cur.execute(queryString)
	seq = cur.fetchone()[1]
	record = SeqRecord(Seq(seq,
                	   IUPAC.unambiguous_dna).translate(),
                   		id=cdsID, name="",
                   		description="")
	records.append(record)

SeqIO.write(records, outFile, "fasta")

listFile.close()
outFile.close()


