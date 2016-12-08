from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from PIL import Image
from StringIO import StringIO

import requests
import json

import collections

from tables import *
from partsdb.tools.Exporters import GenBankExporter
from sqlalchemy import or_

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

def getClassByTablename(tablename):
	for c in Base._decl_class_registry.values():
		if hasattr(c, '__tablename__') and c.__tablename__ == tablename:
			return c
	return None


def getColumnByName(cls, columnName):
	for c in cls.__table__.columns:
		if c.name == columnName:
			return c
	return None



def sortHits(hits, column, nHits):
	sortedHits = sorted( hits, key = lambda x: x[column] )
	return sortedHits[0:nHits+1]

def getGeneHomolog(marpodbSession, geneDBID):
	hits = marpodbSession.query(BlastpHit).\
			filter(BlastpHit.targetID == CDS.id).\
			filter(CDS.id == Gene.cdsID).\
			filter(Gene.dbid == geneDBID ).all()

	if hits:
		hits.sort(key= lambda x: x.eVal)
		return hits[0].proteinName
	else:
		return geneid

def getTopGenes(marpodbSession, StarGene, n):

	geneids = [gene.geneid for gene in StarGene.query.all()]
	topGeneScores = dict(collections.Counter(geneids).most_common(n))

	genes = marpodbSession.query(Gene).filter(Gene.id.in_(topGeneScores.keys()) ).all()

	geneDBIDs = {gene.id : gene.dbid for gene in genes  }

	topGenes = [ ( geneDBIDs[geneid], getGeneHomolog(marpodbSession, geneDBIDs[geneid]), topGeneScores[geneid] ) for geneid in topGeneScores.keys() ]
		
	return sorted(topGenes, key = lambda x: x[2], reverse=True)

def getUserData(StarGene, user, session, marpodbSession):
	userData={}

	if user.is_authenticated:
		geneIds = [ star.geneid for star in StarGene.query.filter(StarGene.userid == user.id)]
	else:
		if not "stars" in session:
			session["stars"] = ""
		geneIds = [i for i in session["stars"].split(':') if i]
		
	userData['starGenes'] = []

	for geneid in geneIds:
		geneDBID = marpodbSession.query(Gene.dbid).filter(Gene.id == geneid).first()
		if geneDBID:
			homolog = getGeneHomolog(marpodbSession, geneDBID)
			userData['starGenes'].append( (geneDBID, homolog ) )

	return userData

def findDataIn(marpodbSession, level, table, queryColumns, equal, returnColumns):

	cls 	= getClassByTablename(table)
	tCls  = cls.__targetclass__

	queryColumns = [ getColumnByName(cls, name) for name in queryColumns ]

	query =  marpodbSession.query(tCls.id, tCls.dbid)

	for column in returnColumns:
		query = query.add_columns(column)

	targets = query.\
				filter(tCls.id == cls.targetID).\
				filter( or_( qC.like(equal) for qC in queryColumns ) ).all()
	if targets:
		return [  (levelID, levelDBID, cols) for levelID, levelDBID, cols in zip( [t[0] for t in targets], [t[1] for t in targets], [ {rc.name: x for rc, x in zip(returnColumns, t[2:]) } for t in targets ] )  ]
	else:
		return {}


def processQuery(marpodbSession, scope, term, columns, nHits):

	# Constructing display columns for the output table 
	displayTables = set( [x.split('.')[1] for x in scope] )

	displayColumns = set()

	for dt in displayTables:
		displayColumns = displayColumns | set( [column.name for column in columns[dt]] )

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
	loci = {}
	for scLevel in scopeDict:
		for scTable in scopeDict[scLevel]:
			
			scColumns = scopeDict[scLevel][scTable]
			newData = findDataIn(marpodbSession, scLevel, scTable, scColumns, term, columns[scTable])

			for hit in newData:

				partID 		= hit[0]
				partDBID	= hit[1]
				cols 		= hit[2]

				fullCols = [''] * len(displayColumns)
				
				for col in cols:
					fullCols[dcMap[col]] = cols[col]

				
				partColumn = getColumnByName(Gene, scLevel+'ID')

				locusDBID = marpodbSession.query(Locus.dbid).filter(Locus.id == Gene.locusID).\
							filter(partColumn == partID).first()[0]

				if locusDBID:
					if not locusDBID in loci:
						loci[locusDBID] = {"genes": {}}

					geneDBID = marpodbSession.query(Gene.dbid).filter(partColumn == partID).first()[0]
					if geneDBID:
						if not geneDBID in loci[locusDBID]["genes"]:
							loci[locusDBID]["genes"][geneDBID] = {"parts": {} }

						if not partDBID in loci[locusDBID]["genes"][geneDBID]["parts"]:
							loci[locusDBID]["genes"][geneDBID]["parts"][partDBID] = {"hits": [] }

						loci[locusDBID]["genes"][geneDBID]["parts"][partDBID]["hits"].append( [scTable] + fullCols )

	if "eVal" in dcMap:
		sortCol = "eVal"
	else:
		sortCol = "name"


	for locusDBID, locus in loci.iteritems():
				# [	z for w in [ genes[gene]     ["trans"][x]["cds"]  [y]["hits"]          for x in genes[gene]     ["trans"]           for y in       genes[gene]     ["trans"][x]["cds"]                         ] for z in w]
		allHits = [ z for w in [ locus["genes"][x]["parts"][y]["hits"] 	       for x in locus["genes"]           for y in       locus["genes"][x]["parts"]                       ] for z in w]

		sortedHits = sortHits(allHits , dcMap[sortCol]+1, nHits)
		locus["topRow"] = [locusDBID] + sortedHits[0][1:]

		# print locusDBID

		for geneDBID, gene in locus["genes"].iteritems():
			# print "\t{0}".format(geneDBID)
						# [z for w in [ genes[gene]["trans"][trans]["cds"][x]["hits"] for x in  genes[gene]["trans"][trans]["cds"]] for z in w ]
			allHits =     [z for w in [ gene["parts"][x]["hits"] 		   					  for x in  gene["parts"]                             ] for z in w ]
			sortedHits = sortHits(allHits, dcMap[sortCol]+1, nHits)
			gene["topRow"] =  [geneDBID] + sortedHits[0][1:]

			for partDBID, part in gene["parts"].iteritems():
				# print "\t\t{0}, {1}".format(partDBID, len(part["hits"]))

				sortedHits = sortHits( part["hits"], dcMap[sortCol]+1, nHits)
				part["topRow"] = [partDBID] + sortedHits[0][1:]
				part["hits"] = sortedHits

	table = {}
	table["header"] = ["id"] + displayColumns
	table["data"] = []

	rowid = 0

	lociList = sorted( loci.values(), key = lambda locus: locus["topRow"][dcMap[sortCol]+1] )

	for locus in lociList:
		rowid += 1
		locusid = rowid

		row = {"rowid": rowid, "pid": "none", "cols": locus["topRow"], "level" : "locus"}
		table["data"].append(row)

		geneList = sorted( locus["genes"].values(), key = lambda gene: gene["topRow"][dcMap[sortCol]+1] )

		for gene in geneList:
			rowid += 1
			geneid = rowid

			row = {"rowid": rowid, "pid": locusid, "cols": gene["topRow"], "level" : "gene"}
			table["data"].append(row)

			partList = sorted( gene["parts"].values(), key = lambda part: part["topRow"][dcMap[sortCol]+1] )

			for part in partList:
				rowid += 1
				partid = rowid

				row = {"rowid": rowid, "pid": geneid, "cols": part["topRow"], "level" : "part"}
				table["data"].append(row)

				for hit in part["hits"]:
					rowid += 1
					row = {"rowid": rowid, "pid": partid, "cols": hit, "level" : "hit"}
					table["data"].append(row)
	return table


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


def getGeneCoordinates(marpodbSession, locusid):
	
	genes = marpodbSession.query(Gene).filter(Gene.locusID == locusid).all()

	exporter = GenBankExporter(None)


	response = {}
	
	response['genes'] = {}

	for gene in genes:
		record =  exporter.export(gene, None)
		response['genes'][record.id] = {'strand' : gene.locusStrand, 'features' : {}}

		for feature in record.features:
			try:
				locations = feature.location.parts
			except:
				locations = [feature.location]

	
			coordinates = ";".join([ "{0}:{1}".format(part.start, part.end) for part in locations  ])

			response['genes'][record.id]['features'][feature.id] = coordinates
		
		if not 'seq' in response:
			print 'Sequence of ', record.id
			if response['genes'][record.id]['strand'] == 1:
				response['seq'] = record.seq
			else:
				response['seq'] = record.seq.reverse_complement()
	
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

