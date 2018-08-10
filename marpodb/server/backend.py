from ..core import *

import os
import json
import requests
import collections

from PIL import Image
from StringIO import StringIO
from sqlalchemy import or_
from Bio.Alphabet import IUPAC

marpodb = MarpoDB("/marpodb4")

def getGeneHomolog(marpodbSession, cdsDBID):
	hit = marpodbSession.query(BlastpHit).\
			filter(BlastpHit.targetID == CDS.id).\
			filter(CDS.dbid == cdsDBID).order_by( BlastpHit.eVal ).first()

	if hit:
		return hit.description
	else:
		return cdsDBID

def findDataIn(marpodbSession, table, queryColumns, equal, returnColumns, dataset):
	cls 			= getClassByTablename(table)
	queryColumns	= [ getColumnByName(cls, name) for name in queryColumns ]
	returnColumns	= [ getColumnByName(cls, columnName) for columnName in returnColumns ] 

	if issubclass( cls, AnnotationMixIn ):
		tCls    = cls.__targetclass__
		
		partIDColumn = getColumnByName(Gene, "{0}ID".format(tCls.__tablename__))
		
		query	= marpodbSession.query(tCls.dbid, Gene.dbid, *returnColumns).filter(Gene.datasetID == Dataset.id).filter(Dataset.name.in_(dataset))
		hits	= query.filter(partIDColumn == tCls.id).filter(tCls.id == cls.targetID).\
						filter( or_( qC.ilike("{0}%".format(equal)) for qC in queryColumns ) ).all()
		if hits:
			return [ ( hit[0], hit[1], {rc.name: val for rc, val in zip(returnColumns, hit[2:])} ) for hit in hits ]

	if issubclass( cls, PartMixIn ):
		partIDColumn = getColumnByName(Gene, "{0}ID".format(cls.__tablename__))
		query	= marpodbSession.query( cls.dbid, Gene.dbid ).filter(Gene.datasetID == Dataset.id).filter(Dataset.name.in_(dataset))
		hits	= query.filter(partIDColumn == cls.id).filter( or_( qC.ilike("{0}%".format(equal)) for qC in queryColumns ) ).all()
		if hits:
			return [ (hit[0], hit[1], {}) for hit in hits ]

	return []

def processQuery(marpodbSession, scope, term, dataset, displayColumns, nHits):

	# Creating a scope map
	scopeDict = {}

	for item in scope:
		tabs = item.split('.')
		scTable = tabs[0]
		scColumn = tabs[1]

		if not scTable in scopeDict:
			scopeDict[scTable] = []
		scopeDict[scTable].append(scColumn)

	# Processing queries
	genes = {}
	for scTable in scopeDict:
		scColumns = scopeDict[scTable]
		newData = findDataIn(marpodbSession, scTable, scColumns, term, displayColumns, dataset)
		for hit in newData:
			partDBID = hit[0]
			geneDBID = hit[1]
			cols = { col: '' for col in displayColumns }
			cols.update(hit[2])
			cols["dbid"] = partDBID


			if not geneDBID in genes:
				genes[geneDBID] = []

			genes[geneDBID].append( cols )

	for geneDBID in genes:
		genes[geneDBID].sort(key = lambda hit: hit["eVal"])

	sortedKeys = sorted( genes, key = lambda k: genes[k][0]["eVal"] )

	table = {}
	table["header"] = ["dbid"] + displayColumns
	table["data"] = []

	rowid = 0
	for geneDBID in sortedKeys:
		rowid += 1
		geneid = rowid
		
		cols = genes[geneDBID][0].copy()
		cols["dbid"] = geneDBID

		row = {"rowid": rowid, "pid": "none", "cols": cols }
		table["data"].append(row)

		for hit in genes[geneDBID][0:5]:
			rowid +=1
			row = {"rowid": rowid, "pid": geneid, "cols": hit }
			table["data"].append(row)

	for row in table["data"]:
		print row
	return table

def parseInterProHTML(html):
	mainContentString = '<div class="grid_19 omega main-content">'
	bottom = mainContentString.join( html.split(mainContentString)[1:] )

	nobody = bottom.split('</body>')[0]

	html = '<div class="container_24">\n\
			<div class="grid_24 clearfix" id="content" >\n\
			<div class="grid_19 omega main-content">' + nobody

	return html.replace('resources/', '/static/interpro/resources/')

def getBlastpHits(marpodbSession, cdsDBID):
	returnTable = {'rows' : [], 'maxLen' : -1}
	hits = marpodbSession.query(BlastpHit).filter(BlastpHit.targetID==CDS.id).filter(CDS.dbid == cdsDBID).all()

	if hits:
		for hit in hits:
			coordinates = [ [hit.start, hit.end, hit.tStart, hit.tEnd, hit.identity] ]
			returnTable["maxLen"] = max(returnTable["maxLen"], hit.tLen)
			
			row = {}
			row['uniID'] = hit.refID
			row['tLen'] = hit.tLen
			row['qLen'] = hit.qLen
			row['proteinName'] = hit.description
			row['origin'] = hit.origin
			row['eVal'] = hit.eVal
			row['coordinates'] = coordinates

			returnTable['rows'].append(row)

		returnTable['maxLen'] = max(returnTable["maxLen"], hits[0].qLen)
	
	returnTable['rows'].sort(key =lambda x: x['eVal'])

	return returnTable

def getDbxRef(marpodbSession, cdsDBID):
	hits = marpodbSession.query(DbxRef.origin, DbxRef.description).filter(DbxRef.targetID==CDS.id).filter(CDS.dbid == cdsDBID).all()
	
	links = {
		"Phytozome MP-Tak 3.1" : "https://phytozome.jgi.doe.gov/pz/portal.html#!results?search=0&crown=1&star=1&method=5040&searchText={0}&offset=0",
		"MarpolBase" : "http://marchantia.info/nomenclature/{0}"
	}

	if hits:
		hits = { hit[0]: [ hit[1], links[hit[0]].format(hit[1]) ] for hit in hits }

	return hits

def getCDSDetails(marpodbSession, cdsDBID):

	response = {}
	response['blastp'] = getBlastpHits(marpodbSession, cdsDBID)
	response['dbxref'] = getDbxRef(marpodbSession, cdsDBID)
	return response

def getGeneDetails(marpodbSession, dbid):
	gene = marpodbSession.query(Gene).filter(Gene.dbid == dbid).first()
	
	if not gene:
		return None

	record = gene.record

	response = { "parts": {}, "seq": record.seq }

	for feature in record.features:
		try:
			locations = feature.location.parts
		except:
			locations = [feature.location]

		coordinates = ";".join([ "{0}:{1}".format(part.start, part.end) for part in locations  ])
		
		response["parts"][feature.id] = coordinates

	return response

def splitBLASTString(s,n):
	return [s[ start:start+n ].replace(" ", "&nbsp") for start in range(0, len(s), n) ]

def parseBlastResult(data, lineLenght = 60):
	blastResult = json.loads(data)
	rows = []

	for hit in blastResult["BlastOutput2"][0]["report"]["results"]["search"]["hits"]:
		row = {}
		title = hit["description"][0]["title"]
		row["dbid"] = title.split()[0]

		row["qLen"]			= blastResult["BlastOutput2"][0]["report"]["results"]["search"]["query_len"]
		row["tLen"]			= hit["len"]
		row["coordinates"]	= []
		row["qseq"]			= []
		row["hseq"]			= []
		row["midline"]		= []
		row["evalue"]		= min( hsp["evalue"] for hsp in hit["hsps"] )
		row["bit_score"]	= min( hsp["bit_score"] for hsp in hit["hsps"] )

		for hsp in hit["hsps"]:
			identity = "{0:.2f}".format(100*float(hsp["identity"]) / row["qLen"])
			coverage = "{0:.2f}".format( float(hsp["align_len"]-hsp["gaps"]) / row["qLen"])
			
			if blastResult["BlastOutput2"][0]["report"]["program"] == "blastx":
				identity = "{0:.2f}".format( 3 * float(identity) )
				coverage = "{0:.2f}".format( 3 * float(coverage) )

			row["qseq"].append( hsp["qseq"] )
			row["hseq"].append( hsp["hseq"] )
			row["midline"].append( hsp["midline"] )
			row["coordinates"].append( [hsp["query_from"], hsp["query_to"], hsp["hit_from"], hsp["hit_to"], identity] )

		row["qseq"] = splitBLASTString( "...".join(row["qseq"]), lineLenght )
		row["hseq"] = splitBLASTString( "...".join(row["hseq"]), lineLenght )
		row["midline"] = splitBLASTString( "...".join(row["midline"]), lineLenght )

		row["len"] = max(row["qLen"] + row["coordinates"][0][2], row["tLen"])

		rows.append(row)
	return rows

def parseHMMResult(tableFileName):
	tableFile = open(tableFileName)
	
	if not tableFile:
		print 'NO FILE', tableFileName

	rows = []

	for line in tableFile:
		if not line.startswith('#'):
			tabs = line.split()
			row = {}
			row['dbid'] = tabs[0]
			row['eVal'] = tabs[4]
			row['score'] = tabs[5]
			row['bias'] = tabs[6]
			rows.append(row)

	return rows

def recodeExport(geneName, seq, seqType):
	record = SeqRecord( Seq( seq, IUPAC.unambiguous_dna ), name = geneName, description = "Generated by MarpoDB", id="M. polymorpha")
	
	strand = 1;
	location = FeatureLocation(11, len(seq)-11, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type=seqType, id="Part" ) )
	
	location = FeatureLocation(0, 6, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="BsaI_F", id="BsaI") )
	
	location = FeatureLocation(7, 11, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="misc_recomb", id="Front" ) )

	strand = -1;
	location = FeatureLocation(len(seq)-6, len(seq), strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="BsaI_R", id="BsaI" ) )
	

	location = FeatureLocation(len(seq)-11, len(seq)-7, strand = strand);
	record.features.append( SeqFeature( location = location, strand = strand, type="misc_recomb", id="Front" ) )

	return record

def recfind(pattern, string, where_should_I_start=0):
	# Save the result in a variable to avoid doing the same thing twice
	pos = string.find(pattern, where_should_I_start)
	if pos == -1:
		# Not found!
		return []
	# No need for an else statement
	return [pos] + recfind(pattern, string, pos + len(pattern))


def generateNewMap(User, staticPath):
	users = User.query.all()
	places = [ user.affiliation for user in users ]
	markers = []

	for place in places:
		print "\n" + place 
		data={"query" : place, "key" : "AIzaSyC3bpj_TeBv6EJgf3zCRd7I__RjKL-hTyU"}
		r = requests.get("https://maps.googleapis.com/maps/api/place/textsearch/json", params = data)
		print r.url
		print r.json
		if r.status_code == requests.codes.ok and r.json()["results"]:
			markers.append(u"{0},{1}".format(r.json()["results"][0]["geometry"]["location"]["lat"], r.json()["results"][0]["geometry"]["location"]["lng"]))	

	with open( os.path.join(staticPath, "json/map_style.json") ) as dataFile:    
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
		mapImage.save(os.path.join(staticPath, "img/map.png"),"PNG")

def getTopGenes(marpodbSession, StarGene, n):

	cdsids = [star.cdsdbid for star in StarGene.query.all()]
	topCDSScores = dict(collections.Counter(cdsids).most_common(n))

	topCDS = [ (cdsdbid, getGeneHomolog(marpodbSession, cdsdbid), topCDSScores[cdsdbid] ) for cdsdbid in topCDSScores.keys() ]
		
	return sorted(topCDS, key = lambda x: x[2], reverse=True)