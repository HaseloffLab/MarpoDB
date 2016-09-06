from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from PIL import Image
from StringIO import StringIO

import requests
import json

import collections

def splitString(s,n):
	return [s[ start:start+n ] for start in range(0, len(s), n) ]

def parseBlastResult(data, lineLenght = 60):
	print "BLAST:"
	print data
	blastResult = json.loads(data)
	rows = []
	for hit in blastResult["BlastOutput2"]["report"]["results"]["search"]["hits"]:
		row = {}
		title = hit["description"][0]["title"]
		row["cdsName"] = '.'.join( title.split()[0].split('_c0_'))
	
		row.update(hit["hsps"][0])

		print "identity: ", row["identity"]
		print "align_len: ", row["align_len"]

		row["identity"] = "{0:.2f}".format(float(row["identity"]) / blastResult["BlastOutput2"]["report"]["results"]["search"]["query_len"])
		row["coverage"] = "{0:.2f}".format( float(row["align_len"]-row["gaps"]) / blastResult["BlastOutput2"]["report"]["results"]["search"]["query_len"])

		row["qseq"] = splitString(row["qseq"], lineLenght )
		row["hseq"] = splitString(row["hseq"], lineLenght )
		row["midline"] = [ s.replace(" ", "&nbsp") for s in splitString(row["midline"], lineLenght )]

		rows.append(row)
	return rows

def generateNewMap(User):

	users = User.query.all()
	places = [ user.affiliation for user in users ]
	markers = []

	print "PLACES:"	

	for place in places:
		print "\n" + place 
		data={"query":place, "key":"AIzaSyC3bpj_TeBv6EJgf3zCRd7I__RjKL-hTyU"}
		r = requests.get("https://maps.googleapis.com/maps/api/place/textsearch/json", params = data)
		print r.url
		print r.json
		if r.status_code == requests.codes.ok:
			markers.append(u"{0},{1}".format(r.json()["results"][0]["geometry"]["location"]["lat"], r.json()["results"][0]["geometry"]["location"]["lng"]))	

	with open('server/static/json/map_style.json') as dataFile:    
	    data = json.load(dataFile)
	
	print markers
	
	data["style"] = [ ("|".join( key + ':' + str(style[key]) for key in sorted(style.keys()) )) for style in data["style"] ]
	data["markers"] = [ u"size:tiny|color:green|{}".format(marker) for marker in markers ]

	r = requests.get("http://maps.googleapis.com/maps/api/staticmap", params = data, stream=True)
	print r.url
	if r.status_code == requests.codes.ok:
		mapImage = Image.open(StringIO(r.content))
		box = (0, 110, 1020, 775)
		mapImage = mapImage.crop(box)
		mapImage.save("server/static/img/map.png","PNG")


def compositeName(name):

	cdsName = ""
	transName = ""

	barsplit = name.split('|')
	if len(barsplit) == 2:
		cdsName = barsplit[1]

	dotSplit = barsplit[0].split('.')
	if len(dotSplit) == 2:
		transName = dotSplit[1]

	return [dotSplit[0], transName, cdsName]

def findDataIn(cur, level, table, queryColumns, equal, returnColumns):
	sReturnColumns = ','.join([ table + '.' + col for col in returnColumns ])
	sQueryColumns = " OR".join( [ " {0}.{1} ILIKE \'%{2}%\'".format(table, column, equal) for column in queryColumns ] )
	queryString = "SELECT {0}.name, {1} FROM {0} JOIN {2} ON ({0}.id = {2}.{0}_id) WHERE {3}".format(level, sReturnColumns, table, sQueryColumns)
	cur.execute(queryString)
	rows = cur.fetchall()
	if rows:
		return [  (levelName, cols) for levelName, cols in zip( [x[0] for x in rows], [ {rc: x for rc, x in zip(returnColumns, row[1:]) } for row in rows ] )  ]
	else:
		return {}



def sortHits(hits, column, nHits):
	sortedHits = sorted( hits, key = lambda x: x[column] )
	return sortedHits[0:nHits+1]

def getGeneHomolog(cur, geneid):
	queryString = "SELECT protein_name, e_val FROM blastp_hit JOIN cds ON blastp_hit.cds_id = cds.id WHERE cds.name LIKE '{0}.%'".format(geneid)
	cur.execute(queryString)
	hits = cur.fetchall()

	if hits:
		hits.sort(key= lambda x: x[1])
		return hits[0][0]
	else:
		return geneid

def getTopGenes(cur, StarGene, n):

	geneids = [gene.geneid for gene in StarGene.query.all()]
	topGeneScores = dict(collections.Counter(geneids).most_common(n))

	print "GENEIDS: ", geneids

	cur.execute("SELECT id, name FROM gene WHERE id IN %s", ( tuple(topGeneScores.keys()),) )
	geneNames = dict(cur.fetchall())

	topGenes = [ ( geneNames[geneid], getGeneHomolog(cur, geneNames[geneid]), topGeneScores[geneid] ) for geneid in topGeneScores.keys() ]


	# for geneId in topGeneIds:
	# 	print geneId
	# 	cur.execute("SELECT name FROM gene WHERE id=%s", (geneId[0], ))
	# 	geneName = cur.fetchone()[0]
	# 	topGenes.append( ( geneName, getGeneHomolog(cur, geneName), geneId[1]) )
		
	return sorted(topGenes, key = lambda x: x[2], reverse=True)

def getUserData(cur, StarGene, user, session):
	userData={}

	if user.is_authenticated:
		geneIds = [ star.geneid for star in StarGene.query.filter(StarGene.userid == user.id)]
	else:
		if not "stars" in session:
			session["stars"] = ""
		geneIds = [i for i in session["stars"].split(':') if i]
		
	userData['starGenes'] = []

	for geneid in geneIds:
		cur.execute("SELECT name FROM gene WHERE id = %s", (geneid,) )
		geneName = cur.fetchone()[0]
		userData['starGenes'].append( (geneName, getGeneHomolog(cur, geneName) ) )

	return userData

def processQuery(cur, scope, term, columns, nHits):

	# Constructing display columns for the output table 
	displayTables = set( [x.split('.')[1] for x in scope] )

	displayColumns = set()
	
	for dt in displayTables:
		displayColumns = displayColumns | set(columns[dt])
		print dt, columns[dt], displayColumns

	displayColumns = list(displayColumns)
	
	dcMap = { v:k for k, v in enumerate(displayColumns) }
	
	# Creating a scope map
	scopeDict = {}

	for item in scope:
		scLevel = item.split('.')[0]
		scTable = item.split('.')[1]
		scColumn = item.split('.')[2]

		if not scLevel in scopeDict:
			scopeDict[scLevel] = {}
			scopeDict[scLevel][scTable] = []
			scopeDict[scLevel][scTable].append(scColumn)
		elif not scTable in scopeDict[scLevel]:
			scopeDict[scLevel][scTable] = []
			scopeDict[scLevel][scTable].append(scColumn)
		else:
			scopeDict[scLevel][scTable].append(scColumn)

	# Processing queries
	genes = {}
	for scLevel in scopeDict:
		for scTable in scopeDict[scLevel]:
			
			queryColumns = scopeDict[scLevel][scTable]
			newData = findDataIn(cur, scLevel, scTable, queryColumns, term, columns[scTable])

			for hit in newData:
				name = compositeName(hit[0])
				cols = hit[1]
				fullCols = [''] * len(displayColumns)

				for col in cols:
					fullCols[dcMap[col]] = cols[col]

				geneName = name[0]

				if not geneName in genes:
					genes[geneName] = { "hits": [], "trans": {} }
				
				if name[1] == "":
					genes[geneName]["hits"].append( [scTable] + fullCols)
				
				else:
					transName = geneName+ '.'+ name[1]
					if not transName in genes[geneName]["trans"]:
						genes[geneName]["trans"][transName] = {"hits": [], "cds" : {} }
					if name[2] == "":
						genes[geneName]["trans"][transName]["hits"].append([scTable] + fullCols)
					else:
						cdsName = transName + '|' + name[2]
						if not cdsName in genes[geneName]["trans"][transName]["cds"]:
							genes[geneName]["trans"][transName]["cds"][cdsName] = {"hits": []}
						genes[geneName]["trans"][transName]["cds"][cdsName]["hits"].append([scTable] + fullCols)

	if "e_val" in dcMap:
		sortCol = "e_val"
	else:
		sortCol = "name"

	# Sorting the table

	for gene in genes:
		allHits = genes[gene]["hits"] + [z for w in [ genes[gene]["trans"][x]["hits"] for x in genes[gene]["trans"] ] for z in w] + [z for w in [ genes[gene]["trans"][x]["cds"][y]["hits"] for x in genes[gene]["trans"] for y in genes[gene]["trans"][x]["cds"] ] for z in w]
		sortedHits = sortHits(allHits , dcMap[sortCol]+1, nHits)
		genes[gene]["topRow"] = [gene] + sortedHits[0][1:]


		genes[gene]["hits"] = sortHits(genes[gene]["hits"], dcMap[sortCol]+1, nHits)

		for trans in genes[gene]["trans"]:
			allHits = genes[gene]["trans"][trans]["hits"] + [z for w in [ genes[gene]["trans"][trans]["cds"][x]["hits"] for x in  genes[gene]["trans"][trans]["cds"]] for z in w ]
			sortedHits = sortHits(allHits, dcMap[sortCol]+1, nHits)
			genes[gene]["trans"][trans]["topRow"] = [trans] + sortedHits[0][1:]

			genes[gene]["trans"][trans]["hits"] = sortHits(genes[gene]["trans"][trans]["hits"], dcMap[sortCol]+1, nHits)

			for cds in genes[gene]["trans"][trans]["cds"]:
				sortedHits = sortHits(genes[gene]["trans"][trans]["cds"][cds]["hits"], dcMap[sortCol]+1, nHits)
				genes[gene]["trans"][trans]["cds"][cds]["topRow"] = [cds] + sortedHits[0][1:]
				genes[gene]["trans"][trans]["cds"][cds]["hits"] = sortedHits

	geneList = sorted(genes.values(), key = lambda gene: gene["topRow"][dcMap[sortCol]+1])
	
	# Generating output table
	table = {}
	table['header'] = ['id'] + displayColumns
	table['data'] = []
	
	rowid = 0

	for gene in geneList:
		rowid += 1
		geneid = rowid

		row = {"rowid": rowid, "pid": "none", "cols": gene["topRow"], "level" : "gene"}
		table["data"].append(row)

		for hit in gene["hits"]:
			rowid += 1
			row = {"rowid": rowid, "pid": geneid, "cols": hit, "level" : "hit"}
			table["data"].append(row)

		for trans in gene["trans"]:
			rowid += 1
			transid = rowid

			row = {"rowid": rowid, "pid": geneid, "cols": gene["trans"][trans]["topRow"], "level" : "trans"}
			table["data"].append(row)

			for hit in gene["trans"][trans]["hits"]:
				rowid += 1
				row = {"rowid": rowid, "pid": transid, "cols": hit, "level" : "hit"}
				table["data"].append(row)

			for cds in gene["trans"][trans]["cds"]:
				rowid += 1
				cdsid = rowid
				row = {"rowid": rowid, "pid": transid, "cols": gene["trans"][trans]["cds"][cds]["topRow"], "level" : "cds"}
				table["data"].append(row)

				for hit in gene["trans"][trans]["cds"][cds]["hits"]:
					rowid += 1
					row = {"rowid": rowid, "pid": cdsid, "cols": hit, "level" : "hit"}
					table["data"].append(row)


	# for gene in genes:
	# 	rowid = rowid + 1
	# 	geneid = rowid

	# 	row = {'rowid': rowid, 'cols': []}
		
	# 	allHits = genes[gene]["hits"] + [z for w in [ genes[gene]["trans"][x]["hits"] for x in genes[gene]["trans"] ] for z in w] + [z for w in [ genes[gene]["trans"][x]["cds"][y]["hits"] for x in genes[gene]["trans"] for y in genes[gene]["trans"][x]["cds"] ] for z in w]
		
	# 	sortedHits = sortHits(allHits , dcMap[sortCol]+1, nHits)

	# 	row['cols'] = [gene] + sortedHits[0][1:]
	# 	table['data'].append(row)

	# 	for hit in sortHits(genes[gene]["hits"], dcMap[sortCol]+1, nHits):
	# 		rowid = rowid + 1
	# 		row = {'rowid': rowid, 'pid': geneid, 'cols': []}
	# 		row['cols'] =  hit
	# 		table['data'].append(row)

	# 	for trans in genes[gene]["trans"]:
	# 		rowid = rowid + 1
	# 		transid = rowid

	# 		row = {'rowid': rowid, 'pid': geneid, 'cols': []}

	# 		allHits = genes[gene]["trans"][trans]["hits"] + [z for w in [ genes[gene]["trans"][trans]["cds"][x]["hits"] for x in  genes[gene]["trans"][trans]["cds"]] for z in w ]
	# 		sortedHits = sortHits(allHits, dcMap[sortCol]+1, nHits)

	# 		row['cols'] = [trans] + sortedHits[0][1:]
	# 		table['data'].append(row)

	# 		for hit in sortHits(genes[gene]["trans"][trans]["hits"], dcMap[sortCol]+1, nHits):
	# 			rowid = rowid + 1
	# 			row = {'rowid': rowid, 'pid': transid, 'cols': []}
	# 			row['cols'] =  hit
	# 			table['data'].append(row)

	# 		for cds in genes[gene]["trans"][trans]["cds"]:
	# 			rowid = rowid + 1
	# 			cdsid = rowid
	# 			row = {'rowid': rowid, 'pid': transid, 'cols': []}

	# 			sortedHits = sortHits(genes[gene]["trans"][trans]["cds"][cds]["hits"], dcMap[sortCol]+1, nHits)

	# 			row['cols'] = [cds] + sortedHits[0][1:]
	# 			table['data'].append(row)

	# 			for hit in sortedHits:
	# 				rowid = rowid + 1
	# 				row = {'rowid': rowid, 'pid': cdsid, 'cols': []}
	# 				row['cols'] =  hit
	# 				table['data'].append(row)

	return table


def getGeneCoordinates(cur, geneName):
	
	# Retrieve transcript coordinates
	queryString = "SELECT transcript.coordinates, transcript.name FROM transcript WHERE name LIKE \'{0}.seq%\'".format(geneName)
	cur.execute(queryString)
	mrnaCoordinates = cur.fetchall()

	geneStart = 1000000
	geneStop = 0

	response = {}
	response['mrnas'] = {}
	response['cdss'] = {}
	response['gene'] = ()

	if mrnaCoordinates:
		for coordinates in mrnaCoordinates:
			exonStrings = coordinates[0].split(';')
			exons = []
			
			for exonString in exonStrings:
				exon = exonString.split(',')
				start = int(exon[0])
				stop = int(exon[1])
				
				if start < geneStart:
					geneStart = start
				if stop > geneStop:
					geneStop = stop
				exons.append( (start, stop) )
			response['mrnas'][coordinates[1]] = sorted(exons, key = lambda x: x[0])
	else:
		return response

	response['gene'] = (geneStart, geneStop)

	# Retrieve cds coordinates
	queryString = "SELECT cds.coordinates, cds.name FROM cds WHERE name LIKE \'{0}.seq%\'".format(geneName)
	cur.execute(queryString)

	cdsCoordinates = cur.fetchall()

	if cdsCoordinates:
		for coordinates in cdsCoordinates:
			cdsStrings = coordinates[0].split(';')
			cdsParts = []

			for cdsString in cdsStrings:
				cds = cdsString.split(',')
				start = int(cds[0])
				stop = int(cds[1])
				direction = cds[2]

				cdsParts.append( (start, stop, direction) )
			response['cdss'][coordinates[1]] = sorted(cdsParts, key = lambda x: x[0])

	return response


def retrieveHits(cur, hitTable, refTable, name):
	returnTable = {'windowSize':0, 'target':[], 'rows':{}}

	queryString = "SELECT {0}.coordinates, {0}.qLen, {0}.tLen, {0}.protein_name, {0}.gene_name, {0}.origin, {0}.e_val FROM {0} JOIN {1} ON ({0}.{1}_id = {1}.id) WHERE {1}.name = \'{2}\'".format(hitTable, refTable, name)
	
	print queryString

	cur.execute(queryString)
	hits = cur.fetchall()
	
	if hits:
		coordinates = [ [ int(x[0].split(',')[0].split(':')[0]), int(x[0].split(',')[0].split(':')[1]), int(x[0].split(',')[1].split(':')[0]), int(x[0].split(',')[1].split(':')[1]) ] for x in hits ]
			
		xMax = max( [c[0] - c[2] for c in coordinates] )
		xMax = max( [xMax, 0] )

		for x,c in zip (hits, coordinates):
			print (x[1] - c[1]) - ( x[2] - c[3] ), x[1], c[1], x[2], c[3]

		yMax = max( [(x[1] - c[1]) - ( x[2] - c[2] - (c[1] - c[0]) ) for x,c in zip (hits, coordinates)] )
		yMax = max( [yMax, 0] )

		print xMax, yMax

		returnTable['windowSize'] = xMax + hits[0][2] + yMax 
		returnTable['target'] = [xMax, hits[0][2]]

		returnTable['rows'] = [ {'name' : x[3], 'origin': x[5], 'e_val': x[6], 'query': [ c[2] - c[0] + xMax, x[1] ], 'hit' : [ c[0], c[1] ]} for x,c in zip(hits, coordinates) ]

	returnTable['rows'] = sorted(returnTable['rows'], key = lambda x: float( x['e_val'] ))

	return returnTable

def getCDSDetails(cur, cdsName):

	transName = cdsName.split('|')[0]
	response = {}

	response['blastm'] = retrieveHits(cur, 'blastm_hit', 'transcript', transName)
	response['blastp'] = retrieveHits(cur, 'blastp_hit', 'cds', cdsName)

	return response

def exportGB(cur, cdsName):
	try:
		transcriptName = cdsName.split('|')[0]
		geneName = transcriptName.split('.')[0]
	except:
		return None

	# Retrieve gene sequence
	cur.execute("SELECT seq FROM gene WHERE name=%s", (geneName, ))
	geneSeq = cur.fetchone()

	if not geneSeq:
		return None


	record = SeqRecord( Seq( geneSeq[0], IUPAC.unambiguous_dna ), name = geneName, description = "Generated by MarpoDB", id=cdsName )

	cur.execute("SELECT coordinates FROM cds WHERE name=%s", (cdsName,))
	cdsCoordinates = cur.fetchone()

	print cdsCoordinates

	if not( cdsCoordinates and cdsCoordinates[0] ):
		return None

	cdsLocations = []

	cdsString = cdsCoordinates[0].split(';')
	print cdsString
	for part in cdsString:
		tabs = part.split(',')
		start = int(tabs[0])-1
		end = 	int(tabs[1])
		strand = 1 if tabs[2]=='+' else -1

		print start, end, strand

		cdsLocations.append( FeatureLocation(start, end, strand = strand) )
	
	if len(cdsLocations) > 1:
		location = CompoundLocation(cdsLocations)
	else:
		location = cdsLocations[0]

	record.features.append( SeqFeature( location = location, strand = strand, type='CDS', id=transcriptName ) )

	cur.execute("SELECT coordinates FROM transcript WHERE name=%s", (transcriptName,))
	transcriptCoordinates = cur.fetchone()

	if not transcriptCoordinates:
		return None

	exons = []

	
	exonStrings = transcriptCoordinates[0].split(';')
	for exon in exonStrings:
		tabs = exon.split(',')
		start = int(tabs[0])-1
		end = int(tabs[1])

		exons.append( FeatureLocation(start, end, strand = strand) )

	if len(exons) > 1:
		location = CompoundLocation(exons)
	else:
		location = exons[0]

	record.features.append( SeqFeature( location = location, strand = strand, type = 'mRNA', id=cdsName ) )
	return record

