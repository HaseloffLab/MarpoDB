import sys
import argparse
from marpodb.core import *

parser = argparse.ArgumentParser(description='MarpoDB4 Extract')

parser.add_argument('dataset',		type=str,
	choices=['tak', 'cam'],
	help = 'dataset to extract from' )

parser.add_argument('tableName', 	type=str,
	choices = ['gene', 'cds', 'promoter', 'terminator', 'utr3', 'utr5'],
	help = 'parts to extract')

parser.add_argument('-p',			action = 'store_true',
	default = False,
	help = 'get protein sequence, for use with CDS')

parser.add_argument('outputFile', type=str, help = 'output file')

args = parser.parse_args()

marpodb = MarpoDB("/marpodb4")

table = getClassByTablename(args.tableName)


session = marpodb.Session()

if args.tableName == 'gene':
	query = session.query(table).filter(table.datasetID == Dataset.id).filter(Dataset.name == args.dataset)
else:
	geneColumn = getColumnByName(Gene, "{0}ID".format(args.tableName))
	query = session.query(table).filter(geneColumn == table.id).filter(Gene.datasetID == Dataset.id).filter(Dataset.name == args.dataset)

parts = query.all()

marpodb.export(parts, args.outputFile, pep = args.p)

session.close()