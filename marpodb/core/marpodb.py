from partsdb.partsdb import PartsDB

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from tables import *

class MarpoDB(PartsDB):
	def __init__(self, dbName, clean = False):
		PartsDB.__init__(self, "postgresql://" + dbName, Base = Base, clean = clean)
		self.setup(prefix="mpdb")
		
	def export(self, parts, fileName, format = "fasta", pep = False):
		records = []
		if not isinstance(parts, (list,)):
			parts = [parts]

		if pep:
			for part in parts:
				if isinstance(part, (CDS,) ):
					records.append(part.pepRecord)
		else:
			for part in parts:
				if isinstance(part, (PartMixIn, Gene) ):
					records.append(part.record)

		outputFile = open(fileName, 'w')
		SeqIO.write(records, outputFile, format)



