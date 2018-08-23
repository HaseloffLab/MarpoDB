import sys
import os
import argparse

from marpodb.core import *

parser = argparse.ArgumentParser(description='MarpoDB4 Extract')

parser.add_argument('dataset',		type=str,
	choices=['tak', 'cam'],
	help = 'dataset to extract from' )

parser.add_argument('tableName', 	type=str, nargs='+',
	choices = ['gene', 'cds', 'promoter', 'terminator', 'utr3', 'utr5'],
	help = 'parts to extract')

parser.add_argument('-p',			action = 'store_true',
	default = False,
	help = 'get protein sequence, for use with CDS')

parser.add_argument('outputDir', type=str, help = 'output file')

args = parser.parse_args()

marpodb = MarpoDB("/marpodb4")
session = marpodb.Session()

for name in args.tableName:
	print "Exporting {0}...".format(name)
	table = getClassByTablename(name)

	if name == "gene":
		query = session.query(table).filter(table.datasetID == Dataset.id).filter(Dataset.name == args.dataset)
	else:
		geneColumn = getColumnByName(Gene, "{0}ID".format(name))
		query = session.query(table).filter(geneColumn == table.id).filter(Gene.datasetID == Dataset.id).filter(Dataset.name == args.dataset)

	parts = query.all()
	outputPath = os.path.join(args.outputDir, "{0}_{1}{2}.fa".format(args.dataset, name, "_p" if args.p else ""))
	marpodb.export(parts, outputPath, pep = args.p)

session.close()

# table = getClassByTablename(args.tableName)
 

# if args.tableName == 'gene':
# 	query = session.query(table).filter(table.datasetID == Dataset.id).filter(Dataset.name == args.dataset)
# else:
# 	geneColumn = getColumnByName(Gene, "{0}ID".format(args.tableName))
# 	query = session.query(table).filter(geneColumn == table.id).filter(Gene.datasetID == Dataset.id).filter(Dataset.name == args.dataset)

# parts = query.all()

# marpodb.export(parts, args.outputFile, pep = args.p)

# session.close()