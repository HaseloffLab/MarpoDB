import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def saveSequences(rawRecords, fileName, translate = False, format = "fasta"):

	records = []

	for record in rawRecords:
		if translate:
			record = SeqRecord(Seq(record[-1],
	                	IUPAC.ambiguous_dna).translate(),
	                   	id=record[0], name="",
	                   	description='|'.join(record[1:-1]))
		else:
			record = SeqRecord(Seq(record[-1],
	                	IUPAC.ambiguous_dna),
	                   	id=record[0], name="",
	                   	description='|'.join(record[1:-1]))
		records.append(record)


	outFile  = open(fileName, 'w')
	SeqIO.write(records, outFile, format)
	outFile.close()
