from partsdb.partsdb import PartsDB

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from tables import *

def getColumnByName(cls, columnName):
	if columnName in cls.__table__.columns:
		return cls.__table__.columns[columnName]
	else:
		return None
	for c in cls.__table__.columns:
		if c.name == columnName:
			return c
	return None

def getClassByTablename(tablename):
	for c in Base._decl_class_registry.values():
		if hasattr(c, '__tablename__') and c.__tablename__ == tablename:
			return c
	return None

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



