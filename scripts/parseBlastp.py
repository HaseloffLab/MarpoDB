from __future__ import division
import sys
import psycopg2

inFile = open(sys.argv[1])

dbName = sys.argv[2]

conn = psycopg2.connect("dbname={0}".format(dbName))
cur = conn.cursor()

hits = {}

commonPrefix = "mpdb"

def newID(table):
	cur.execute("SELECT COUNT(*) FROM {0}".format(table), )
	nid = cur.fetchone()[0]
	if nid is not None:
		nid = "{0}.{1}.{2}".format(commonPrefix, table, int(nid) )
		return nid
	else:
		return None

def insertRow(table, valuesDict):
	nid = newID(table)

	if nid:
		valuesDict['id'] = nid
		keys = tuple( valuesDict.keys() )
		values = tuple( [ valuesDict[key] for key in keys ] )
		
		keyString = ','.join(keys)

		exeString = "INSERT INTO {0} ({1}) VALUES %s RETURNING id".format(table, keyString)
		cur.execute(exeString, (values,))
		

		cid = cur.fetchone()[0]

		if cid:
			conn.commit()
			return cid

	return None

for line in inFile:
	tabs = line.rstrip().split('\t')
	
	names = ["cdsID", "uniID", "qlen", "slen", "qstart", "qend", "tstart", "tend", "qcovs", "pident", "evalue", "proteinName", "origin",  "geneName"]

	data = dict( zip(names, tabs) )
	data["uniID"] = data["uniID"].split('|')[2]

	name = "{0}_{1}".format(data["cdsID"], data["uniID"])
	
	if not name in hits:
		hits[name] = {}
		hits[name]["cds_id"]			= data["cdsID"]
		hits[name]["uni_id"]			= data["uniID"]
		hits[name]["coverage"]			= data["qcovs"]
		hits[name]["qlen"] 				= data["qlen"]
		hits[name]["tlen"] 				= data["slen"]
		hits[name]["coordinates"]		= ""
		hits[name]["e_val"] 			= float(data["evalue"])
		hits[name]["protein_name"] 		= data["proteinName"]
		hits[name]["gene_name"] 		= data["geneName"]
		hits[name]["origin"] 			= data["origin"]
		
	
	hits[name]["e_val"] = min( hits[name]["e_val"], float(data["evalue"]) )

	coordinates = "{0}:{1},{2}:{3},{4};".format(data["qstart"], data["qend"], data["tstart"], data["tend"], data["pident"])
	hits[name]["coordinates"] += coordinates

inFile.close()

for hitName, hit in hits.iteritems():
	cdsID = hitName.split('_')[0]
	uniID = hitName.split('_')[1]

	cur.execute("SELECT id from cds WHERE id = %s", (cdsID,))
	cdsID = cur.fetchone()[0]

	if cdsID:
		insertRow('blastp_hit', hit)
	else:
		print "Failed to locate {0}".format(hitName.split('_')[0])

conn.commit()
cur.close()
conn.close()
