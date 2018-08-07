import sys
from Bio import SeqIO

from marpodb.core import *

marpodb = MarpoDB('/marpodb4')

inFile = open( sys.argv[1] )

session = marpodb.Session()

records = []

for line in inFile:
	mpid = line.rstrip()

	promoter = session.query(Promoter)\
				.filter(Promoter.id == Gene.promoterID)\
				.filter(DbxRef.targetID == Gene.cdsID)\
				.filter(DbxRef.refID == mpid)\
				.filter(Gene.datasetID == Dataset.id)\
				.filter(Dataset.name == 'tak').first()
	if promoter:
		record = promoter.record
		record.id = mpid

		if not record.seq:
			print "{0} EMPTY SEQ".format(mpid)

		records.append( record )
	else:
		print "{0} NOT FOUND".format(mpid)

inFile.close()

outFile = open(sys.argv[2], 'w')

SeqIO.write(records, outFile, "fasta")