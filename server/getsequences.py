import sys
from partsdb.partsdb import PartsDB
from tables import *

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from partsdb.tools.Exporters import GenBankExporter
from Bio import SeqIO

exportFeature = sys.argv[2]

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)
session = marpodb.Session()

records = []

exporter = GenBankExporter(marpodb)

if exportFeature == 'Prot':

	for cds in session.query(CDS).all():
		feature = SeqFeature( location = exporter.coordinatesToLocation(cds.coordinates) )
		seq = feature.extract(Seq(cds.seq, generic_dna)).translate()
		
		if not (seq[0] == 'M' and seq.find('*') == len(seq)-1):
			print gene.cds.dbid, seq
		else:
			record = SeqRecord( seq = seq, id = cds.dbid, description='Extracted from '+sys.argv[1] )
			records.append(record)

elif exportFeature == 'CDS':
	for cds in session.query(CDS).all():
		feature = SeqFeature( location = exporter.coordinatesToLocation(cds.coordinates) )
		seq = feature.extract(Seq(cds.seq, generic_dna))

		record = SeqRecord( seq = seq, id = cds.dbid, description='Extracted from '+sys.argv[1] )
		records.append(record)

elif exportFeature == 'Gene':
	for gene in session.query(Gene).all():
		seq = ''
		if gene.promoter:
			seq += gene.promoter.seq
		if gene.utr5:
			seq += gene.utr5.seq
		if gene.utr3:
			seq += gene.utr3.seq
		if gene.terminator:
			seq += gene.terminator.seq
		record = SeqRecord( seq = Seq(seq), id = gene.dbid, description='Extracted from '+sys.argv[1] )
		records.append(record)

else:
	print "Unknown exportFeature: '{0}'".format(exportFeature)
	sys.exit() 


outputFile = open('data/{0}.fa'.format(exportFeature), 'w')

SeqIO.write(records, outputFile, 'fasta')
print "Exported {0} sequences to data/{1}.fa".format(len(records), exportFeature)

session.close()
