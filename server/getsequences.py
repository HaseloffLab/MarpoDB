import sys
from partsdb.partsdb import PartsDB
from tables import *

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from partsdb.tools.Exporters import GenBankExporter
from Bio import SeqIO

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

records = []

exporter = GenBankExporter(marpodb)

for gene in marpodb.session.query(Gene).all():
	feature = SeqFeature( location = exporter.coordinatesToLocation(gene.cds.coordinates) )
	seq = feature.extract(Seq(gene.cds.seq, generic_dna)).translate()
	if not (seq[0] == 'M' and seq.find('*') == len(seq) - 1):
		print gene.cds.dbid, seq
	else:
		record = SeqRecord( seq = seq, id = gene.cds.dbid )
		records.append(record)

outputFile = open('CDSs.fa', 'w')

SeqIO.write(records, outputFile, 'fasta')
