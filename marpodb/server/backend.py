from ..core import *

from sqlalchemy import or_
import pandas as pd


marpodb = MarpoDB("/marpodb4")

def getColumnByName(cls, columnName):
	if columnName in cls.__table__.columns:
		return cls.__table__.columns[columnName]
	else:
		return None
	for c in cls.__table__.columns:
		if c.name == columnName:
			return c
	return None

def getClassByTablename(tablename):
	for c in Base._decl_class_registry.values():
		if hasattr(c, '__tablename__') and c.__tablename__ == tablename:
			return c
	return None


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