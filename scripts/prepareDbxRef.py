from __future__ import division
import sys

headers = ["dbid", "origin", "start", "end", "description", "refID", "eVal"]

inFile = open( sys.argv[1] ) 

annotated = []

print "\t".join(headers)

for line in inFile:
	tabs 	= line.split()
	dbid 	= tabs[0]
	phtz 	= tabs[1]
	qlen 	= int(tabs[2])
	slen 	= int(tabs[3])
	lenght 	= int(tabs[4])
	ident	= float(tabs[5])

	alias = phtz.split('.')[0]
	cov   = 100 * lenght / qlen

	if ident > 95 and cov > 95 and qlen == slen and not dbid in annotated:
		annotated.append(dbid)
		print "\t".join([dbid, "Phytozome MP-Tak 3.1", "0", str(qlen), alias, alias, "0"])
