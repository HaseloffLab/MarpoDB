import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def saveSequences(rawRecords, fileName, genomic = False, translate = False, format = "fasta"):

	records = []
	for record in rawRecords:
		
		id = record.dbid
		
		try:
			coordinates = record.exons
		except:
			coordinates = None

		sequence = record.seq

		if coordinates and not genomic:
			exons = [ (ex.split(',')[0], ex.split(',')[1]) for ex in coordinates.split('+') ]
			seq = ''
			for exon in exons:
				seq += sequence[ exon[0]-1 : exon[1] ]
			sequence = seq

		if translate:
			record = SeqRecord(Seq(sequence,
	                	IUPAC.ambiguous_dna).translate(),
	                   	id=id, name="",
	                   	description="")
		else:
			record = SeqRecord(Seq(sequence,
	                	IUPAC.ambiguous_dna),
	                   	id=id, name="",
	                   	description="")
		records.append(record)


	outFile  = open(fileName, 'w')
	SeqIO.write(records, outFile, format)
	outFile.close()
