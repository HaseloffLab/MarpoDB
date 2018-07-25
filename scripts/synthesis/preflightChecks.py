from utils import maxSynLen, minPromoterLen, minSynLen

from marpodb.client import MarpoDBClient
from marpodb.client.tables import *

client = MarpoDBClient("/marpodbtak")
session = client.Session()

inputList = open("../../data/synthesis/21_05_2018_TF_SYNTHESIS.tsv")

keys = checks.keys()

print "\t".join(["alias", "tff", "info_name"] + keys)

for line in inputList:
	if not line.startswith("#"):
		
		tabs = line.split("\t")

		alias = tabs[0]
		tff = tabs[1]
		info_name = tabs[2].rstrip()

		gene = session.query(Gene).filter(Gene.alias == alias).first()

		row = [ alias, tff, info_name ]

		for key in checks.keys():
			row.append( str( checks[key](gene) ) )

		print "\t".join( row )


