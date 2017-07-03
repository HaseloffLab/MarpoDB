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
from sqlalchemy import or_, desc

def splitString(s,n):
	return [s[ start:start+n ] for start in range(0, len(s), n) ]

def recfind(pattern, string, where_should_I_start=0):
    # Save the result in a variable to avoid doing the same thing twice
    pos = string.find(pattern, where_should_I_start)
    if pos == -1:
        # Not found!
        return []
    # No need for an else statement
    return [pos] + recfind(pattern, string, pos + len(pattern))

def parseHMMResult(tableFileName, session):
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
			try:
				row['locusdbid'] = session.query(Locus.dbid).filter(CDS.dbid == row["dbid"]).filter( Gene.cdsID == CDS.id ).filter(Locus.id == Gene.locusID).first()[0]
			except:
				continue
			rows.append(row)

	return rows


def parseBlastResult(data, session, lineLenght = 60):
	print "BLAST:"
	# print data
	blastResult = json.loads(data)
	rows = []

	for hit in blastResult["BlastOutput2"][0]["report"]["results"]["search"]["hits"]:
		row = {}
		title = hit["description"][0]["title"]
		row["dbid"] = title.split()[0]
		
		# Quick fix...
		try:
			if 'cds' in row["dbid"]:
				row['locusdbid'] = session.query(Locus.dbid).filter(CDS.dbid == row["dbid"]).filter( Gene.cdsID == CDS.id ).filter(Locus.id == Gene.locusID).first()[0]
			elif 'gene' in row["dbid"]:
				row['locusdbid'] = session.query(Locus.dbid).filter(Gene.dbid == row['dbid']).filter(Gene.locusID == Locus.id).first()[0]
			row.update(hit["hsps"][0])
		except:
			continue

		print "identity: ", row["identity"]
		print "align_len: ", row["align_len"]

		row["identity"] = "{0:.2f}".format(100*float(row["identity"]) / blastResult["BlastOutput2"][0]["report"]["results"]["search"]["query_len"])
		row["coverage"] = "{0:.2f}".format( float(row["align_len"]-row["gaps"]) / blastResult["BlastOutput2"][0]["report"]["results"]["search"]["query_len"])
		if blastResult["BlastOutput2"][0]["report"]["program"] == "blastx":
			row["identity"] = "{0:.2f}".format( 3 * float(row["identity"]) )
			row["coverage"] = "{0:.2f}".format( 3 * float(row["coverage"]) )

		row["qseq"] = splitString(row["qseq"], lineLenght )
		row["hseq"] = splitString(row["hseq"], lineLenght )
		row["midline"] = [ s.replace(" ", "&nbsp") for s in splitString(row["midline"], lineLenght )]
		row["coordinates"] = [ [ row["query_from"], row["query_to"], row["hit_from"], row["hit_to"], row["identity"]]  ]
		row["qLen"] = blastResult["BlastOutput2"][0]["report"]["results"]["search"]["query_len"]
		row["tLen"] = hit["len"]
		row["len"] = max(row["qLen"], row["tLen"])

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

def getGeneHomolog(marpodbSession, cdsDBID):
	hit = marpodbSession.query(BlastpHit).\
			filter(BlastpHit.targetID == CDS.id).\
			filter(CDS.dbid == cdsDBID).order_by( BlastpHit.eVal ).first()

	if hit:
		return hit.description
	else:
		return cdsDBID

def getTopGenes(marpodbSession, StarGene, n):

	cdsids = [star.cdsdbid for star in StarGene.query.all()]
	topCDSScores = dict(collections.Counter(cdsids).most_common(n))

	topCDS = [ (cdsdbid, getGeneHomolog(marpodbSession, cdsdbid), topCDSScores[cdsdbid] ) for cdsdbid in topCDSScores.keys() ]
		
	return sorted(topCDS, key = lambda x: x[2], reverse=True)

def getUserData(userDB, StarGene, user, session, marpodbSession):
	userData={}

	if user.is_authenticated:
		stars = StarGene.query.filter(StarGene.userid == user.id).all()
		
		userData['starGenes'] = {}

		for star in stars:
			if not star.genename:
				star.genename = getGeneHomolog(marpodbSession, star.cdsdbid)
				userDB.session.add(star)
			userData['starGenes'][star.cdsdbid] = star.genename
		userDB.session.commit()
	return userData

def findDataIn(marpodbSession, table, queryColumns, equal, returnColumns):
	cls 	= getClassByTablename(table)
	queryColumns  = [ getColumnByName(cls, name) for name in queryColumns ]
	returnColumns = [ getColumnByName(cls, columnName) for columnName in returnColumns ] 

	if issubclass( cls, AnnotationMixIn ):
		tCls    = cls.__targetclass__
		query = marpodbSession.query(tCls)
		for column in returnColumns:
			query = query.add_columns( column )

		targets = query.filter(tCls.id == cls.targetID).\
						filter( or_( qC.ilike(equal+'%') for qC in queryColumns ) ).all()
		print targets
		if targets:
			return [ ( partDBID, geneDBID, cols ) for partDBID, geneDBID, cols in zip(
						[ target[0].dbid for target in targets ],
						[ target[0].gene[0].dbid for target in targets ],
						[ {rc.name: val for rc, val in zip( returnColumns, target[1:])} for target in targets ]
					) ]
	else:
		if issubclass(cls, PartMixIn):
			print "CLS: ", cls.__tablename__
			print "Columns: ", queryColumns
			query = marpodbSession.query(cls.dbid, cls)
			targets = query.filter( or_(qC.ilike(equal+'%') for qC in queryColumns) ).all()
			if targets:
				print "return: ", returnColumns
				return [ (partDBID, geneDBID, {}) for partDBID, geneDBID in zip(
							[ target[0] for target in targets ],
							[ target[1].gene[0].dbid for target in targets ]
						) ]
	
	return {}

def processQuery(marpodbSession, scope, term, displayColumns, nHits):

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
		newData = findDataIn(marpodbSession, scTable, scColumns, term, displayColumns)
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

		for hit in genes[geneDBID]:
			rowid +=1
			row = {"rowid": rowid, "pid": geneid, "cols": hit }
			table["data"].append(row)

	for row in table["data"]:
		print row
	return table


	# loci = {}
	# for scLevel in scopeDict:
	# 	for scTable in scopeDict[scLevel]:

	# 		scColumns = scopeDict[scLevel][scTable]

	# 		print scTable, scColumns

	# 		newData = findDataIn(marpodbSession, scLevel, scTable, scColumns, term, columns[scTable])

	# 		print newData

	# 		for hit in newData:

	# 			partID 		= hit[0]
	# 			partDBID	= hit[1]
	# 			cols 		= hit[2]

	# 			fullCols = [''] * len(displayColumns)
				
	# 			for col in cols:
	# 				fullCols[dcMap[col]] = cols[col]

				
	# 			partColumn = getColumnByName(Gene, scLevel+'ID')

	# 			print partColumn, partID

	# 			# locusDBID = marpodbSession.query(Locus.dbid).filter(Locus.id == Gene.locusID).\
	# 			# 			filter(partColumn == partID).first()[0]

	# 			out = marpodbSession.query(Locus.dbid, Gene.dbid).filter( Locus.id == Gene.locusID ).\
	# 						filter(partColumn == partID).first()
				
	# 			if not out:
	# 				continue

	# 			locusDBID = out[0]
	# 			geneDBID  = out[1]

	# 			if locusDBID:
	# 				if not locusDBID in loci:
	# 					loci[locusDBID] = {"genes": {}}

	# 				# geneDBID = marpodbSession. query(Gene.dbid).filter(partColumn == partID).first()[0]
	# 				if geneDBID:
	# 					if not geneDBID in loci[locusDBID]["genes"]:
	# 						loci[locusDBID]["genes"][geneDBID] = {"parts": {} }

	# 					if not partDBID in loci[locusDBID]["genes"][geneDBID]["parts"]:
	# 						loci[locusDBID]["genes"][geneDBID]["parts"][partDBID] = {"hits": [] }

	# 					loci[locusDBID]["genes"][geneDBID]["parts"][partDBID]["hits"].append( [scTable] + fullCols )

	# if "eVal" in dcMap:
	# 	sortCol = "eVal"
	# else:
	# 	sortCol = "name"


	# for locusDBID, locus in loci.iteritems():
	# 			# [	z for w in [ genes[gene]     ["trans"][x]["cds"]  [y]["hits"]          for x in genes[gene]     ["trans"]           for y in       genes[gene]     ["trans"][x]["cds"]                         ] for z in w]
	# 	allHits = [ z for w in [ locus["genes"][x]["parts"][y]["hits"] 	       for x in locus["genes"]           for y in       locus["genes"][x]["parts"]                       ] for z in w]

	# 	sortedHits = sortHits(allHits , dcMap[sortCol]+1, nHits)
	# 	locus["topRow"] = [locusDBID] + sortedHits[0][1:]

	# 	# print locusDBID

	# 	for geneDBID, gene in locus["genes"].iteritems():
	# 		# print "\t{0}".format(geneDBID)
	# 					# [z for w in [ genes[gene]["trans"][trans]["cds"][x]["hits"] for x in  genes[gene]["trans"][trans]["cds"]] for z in w ]
	# 		allHits =     [z for w in [ gene["parts"][x]["hits"] 		   					  for x in  gene["parts"]                             ] for z in w ]
	# 		sortedHits = sortHits(allHits, dcMap[sortCol]+1, nHits)
	# 		gene["topRow"] =  [geneDBID] + sortedHits[0][1:]

	# 		for partDBID, part in gene["parts"].iteritems():
	# 			# print "\t\t{0}, {1}".format(partDBID, len(part["hits"]))

	# 			sortedHits = sortHits( part["hits"], dcMap[sortCol]+1, nHits)
	# 			part["topRow"] = [partDBID] + sortedHits[0][1:]
	# 			part["hits"] = sortedHits

	# table = {}
	# table["header"] = ["id"] + displayColumns
	# table["data"] = []

	# rowid = 0

	# lociList = sorted( loci.values(), key = lambda locus: locus["topRow"][dcMap[sortCol]+1] )

	# for locus in lociList:
	# 	rowid += 1
	# 	locusid = rowid

	# 	row = {"rowid": rowid, "pid": "none", "cols": locus["topRow"], "level" : "locus"}
	# 	table["data"].append(row)

	# 	geneList = sorted( locus["genes"].values(), key = lambda gene: gene["topRow"][dcMap[sortCol]+1] )

	# 	for gene in geneList:
	# 		rowid += 1
	# 		geneid = rowid

	# 		row = {"rowid": rowid, "pid": locusid, "cols": gene["topRow"], "level" : "gene"}
	# 		table["data"].append(row)

	# 		partList = sorted( gene["parts"].values(), key = lambda part: part["topRow"][dcMap[sortCol]+1] )

	# 		for part in partList:
	# 			rowid += 1
	# 			partid = rowid

	# 			row = {"rowid": rowid, "pid": geneid, "cols": part["topRow"], "level" : "part"}
	# 			table["data"].append(row)

	# 			for hit in part["hits"]:
	# 				rowid += 1
	# 				row = {"rowid": rowid, "pid": partid, "cols": hit, "level" : "hit"}
	# 				table["data"].append(row)
	# return table


	# # Sorting the table
	# for gene in genes:
	# 	allHits = genes[gene]["hits"] + [z for w in [ genes[gene]["trans"][x]["hits"] for x in genes[gene]["trans"] ] for z in w] + [z for w in [ genes[gene]["trans"][x]["cds"][y]["hits"] for x in genes[gene]["trans"] for y in genes[gene]["trans"][x]["cds"] ] for z in w]
	# 	sortedHits = sortHits(allHits , dcMap[sortCol]+1, nHits)
	# 	genes[gene]["topRow"] = [gene] + sortedHits[0][1:]


	# 	genes[gene]["hits"] = sortHits(genes[gene]["hits"], dcMap[sortCol]+1, nHits)

	# 	for trans in genes[gene]["trans"]:
	# 		allHits = genes[gene]["trans"][trans]["hits"] + [z for w in [ genes[gene]["trans"][trans]["cds"][x]["hits"] for x in  genes[gene]["trans"][trans]["cds"]] for z in w ]
	# 		sortedHits = sortHits(allHits, dcMap[sortCol]+1, nHits)
	# 		genes[gene]["trans"][trans]["topRow"] = [trans] + sortedHits[0][1:]

	# 		genes[gene]["trans"][trans]["hits"] = sortHits(genes[gene]["trans"][trans]["hits"], dcMap[sortCol]+1, nHits)

	# 		for cds in genes[gene]["trans"][trans]["cds"]:
	# 			sortedHits = sortHits(genes[gene]["trans"][trans]["cds"][cds]["hits"], dcMap[sortCol]+1, nHits)
	# 			genes[gene]["trans"][trans]["cds"][cds]["topRow"] = [cds] + sortedHits[0][1:]
	# 			genes[gene]["trans"][trans]["cds"][cds]["hits"] = sortedHits

	# geneList = sorted(genes.values(), key = lambda gene: gene["topRow"][dcMap[sortCol]+1])
	
	# # Generating output table
	# table = {}
	# table['header'] = ['id'] + displayColumns
	# table['data'] = []
	
	# rowid = 0

	# for gene in geneList:
	# 	rowid += 1
	# 	geneid = rowid

	# 	row = {"rowid": rowid, "pid": "none", "cols": gene["topRow"], "level" : "gene"}
	# 	table["data"].append(row)

	# 	for hit in gene["hits"]:
	# 		rowid += 1
	# 		row = {"rowid": rowid, "pid": geneid, "cols": hit, "level" : "hit"}
	# 		table["data"].append(row)

	# 	for trans in gene["trans"]:
	# 		rowid += 1
	# 		transid = rowid

	# 		row = {"rowid": rowid, "pid": geneid, "cols": gene["trans"][trans]["topRow"], "level" : "trans"}
	# 		table["data"].append(row)

	# 		for hit in gene["trans"][trans]["hits"]:
	# 			rowid += 1
	# 			row = {"rowid": rowid, "pid": transid, "cols": hit, "level" : "hit"}
	# 			table["data"].append(row)

	# 		for cds in gene["trans"][trans]["cds"]:
	# 			rowid += 1
	# 			cdsid = rowid
	# 			row = {"rowid": rowid, "pid": transid, "cols": gene["trans"][trans]["cds"][cds]["topRow"], "level" : "cds"}
	# 			table["data"].append(row)

	# 			for hit in gene["trans"][trans]["cds"][cds]["hits"]:
	# 				rowid += 1
	# 				row = {"rowid": rowid, "pid": cdsid, "cols": hit, "level" : "hit"}
	# 				table["data"].append(row)

	# return table


def getGeneCoordinates(marpodbSession, dbid):
	gene = marpodbSession.query(Gene).filter(Gene.dbid == dbid).first()
	
	if not gene:
		return None

	exporter = GenBankExporter(None)

	record = exporter.export(gene, None)

	response = { "parts": {}, "seq": record.seq }

	for feature in record.features:
		try:
			locations = feature.location.parts
		except:
			locations = [feature.location]

		coordinates = ";".join([ "{0}:{1}".format(part.start, part.end) for part in locations  ])
		
		response["parts"][feature.id] = coordinates

	return response

# def getGeneCoordinates(marpodbSession, locusid):
	
# 	genes = marpodbSession.query(Gene).filter(Gene.locusID == locusid).all()

# 	exporter = GenBankExporter(None)



# 	response = {}
	
# 	response['genes'] = {}

# 	for gene in genes:
# 		record =  exporter.export(gene, None)
		
# 		# I do not know why but if I change the gene.locusStrand variable to 1 then it makes it right, at least for some of them...
# 		response['genes'][record.id] = {'strand' : gene.locusStrand, 'features' : {}}
# 		#response['genes'][record.id] = {'strand' : 1, 'features' : {}}
		
# 		for feature in record.features:
# 			try:
# 				locations = feature.location.parts
# 			except:
# 				locations = [feature.location]
	
# 			coordinates = ";".join([ "{0}:{1}".format(part.start, part.end) for part in locations  ])
			
# 			response['genes'][record.id]['features'][feature.id] = coordinates
		
# 		if not 'seq' in response:
# 			print 'Sequence of ', record.id
# 			if response['genes'][record.id]['strand'] == 1:
# 				response['seq'] = record.seq
# 			else:
# 				response['seq'] = record.seq.reverse_complement()
	
# 	return response

def getBlastpHits(marpodbSession, cdsDBID):
	returnTable = {'rows' : [], 'maxLen' : -1}
	hits = marpodbSession.query(BlastpHit).filter(BlastpHit.targetID==CDS.id).filter(CDS.dbid == cdsDBID).all()

	if hits:
		for hit in hits:
			# coordinateString = hit.coordinates
			# tabs = coordinateString.split(';')
			# tabs = [coordinateString]
			# coordinates = [  [int(tab.split(',')[0].split(':')[0]), int(tab.split(',')[0].split(':')[1]),\
			# 								int(tab.split(',')[1].split(':')[0]), int(tab.split(',')[1].split(':')[1]),\
			# 									float(tab.split(',')[2])] for tab in tabs ]
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
	hit = marpodbSession.query(DbxRef.origin, DbxRef.description).filter(DbxRef.targetID==CDS.id).filter(CDS.dbid == cdsDBID).first()
	hit = { hit[0]: hit[1] } if hit else {}

	return hit

def getCDSDetails(marpodbSession, cdsDBID):

	response = {}
	response['blastp'] = getBlastpHits(marpodbSession, cdsDBID)
	response['dbxref'] = getDbxRef(marpodbSession, cdsDBID)
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

