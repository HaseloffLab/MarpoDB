## Extract parts from the database and domesticate
## Input format:
## [0]name	[1]type [2]batch [3]description [4]dbid  [5]Comments [6]author [7]group [8]reference

import sys
import os

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from utils import *
from marpodb.client import MarpoDBClient
from marpodb.client.tables import *

inList = open( sys.argv[1] )
outPath = sys.argv[3]

if not os.path.isdir( outPath ):
	os.mkdir(outPath)

client = MarpoDBClient("/marpodbtak")
session = client.Session()

parts = {}

for line in inList:
	tabs = line.split('\t')

	name = tabs[0]
	type = tabs[1]
	description = tabs[3]
	author = tabs[6]
	group = tabs[7]
	reference = tabs[8]
	dbid = tabs[4]

	gene = session.query(Gene).filter(Gene.dbid == dbid).first()

	if gene:
		if type == 'PROM5':
			promExtract = extractPROM5(gene, name)
			if promExtract:
				if len(promExtract) == 1:
					parts["PROM5_{0}".format(name)] = domesticate(promExtract[0], 'PROM5', name, description, reference, author, group, gene.alias, addBsaI = True)
				elif len(promExtract) == 2:
					parts['PROM_{0}'.format(name)] = domesticate(promExtract[0], 'PROM', name, description, reference, author, group, gene.alias, addBsaI = True)
					parts['5UTR_{0}'.format(name)] = domesticate(promExtract[1], '5UTR', name, description, reference, author, group, gene.alias, addBsaI = True)
		else:
			rec = SeqRecord( seq = extract[type](gene) )
			parts["{0}_{1}".format(type, name)] = domesticate( rec, type, name, description, reference, author, group, "{0}({1})".format(gene.alias, dbid), addBsaI = True )
	else:
		print "Can't find {0} : {1}".format(name, dbid)
		continue

for part in parts:
	if parts[part]:
		outFile = os.path.join(outPath, part + '.gb')
		SeqIO.write(parts[part], outFile, "genbank")


session.close()

	