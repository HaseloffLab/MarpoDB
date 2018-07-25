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

def queryTable(session, table, queryColumns, equal, returnColumns):
	
	if table in marpodb.classes:
		cls = marpodb.classes[table]
	else:
		return {}

	labels = ["partID"] + returnColumns

	queryColumns  = [ getColumnByName(cls, name) for name in queryColumns ]
	returnColumns = [ getColumnByName(cls, columnName) for columnName in returnColumns ] 

	if issubclass( cls, AnnotationMixIn ):
		tCls    = cls.__targetclass__
		query = session.query(cls.targetID)

		for column in returnColumns:
			query = query.add_columns( column )

		hits = query.filter( or_( qC.ilike(equal+'%') for qC in queryColumns ) ).all()

		if hits:
			hitDataFrame = pd.DataFrame.from_records(hits, columns=labels)
			# hitDataFrame["partDBID"] = session.query()
			
			return hitDataFrame

			# Add columns that contain partDBID and geneDB id. Remove transcript evidence blast to make everything quick.

			# data = {
			# 	"pdbid" : [hit.target.dbid for hit in hits],
			# 	"gdbid" : [hit.target.gene.dbid for ]
			# }


		# if hits:

		# 	partDBIDs = [hit.target.dbid for hit in hits]
		# 	geneDBIDs = [hit.target.gene[0].dbid for hit in hits ]
		# 	cols = [{ rc: getattr(hit, rc) for rc in returnColumns } for hit in hits]
 	# 		return [ (pdbid, dbid, col) for pdbid, dbid, col in zip( partDBIDs, geneDBIDs, cols ) ]
			# return [ ( partDBID, geneDBID, cols ) for partDBID, geneDBID, cols in zip(
			# 			[ target[0].dbid for target in targets ],
			# 			[ target[0].gene[0].dbid for target in targets ],
			# 			[ {rc.name: val for rc, val in zip( returnColumns, target[1:])} for target in targets ]
			# 		) ]
	else:
		if issubclass(cls, PartMixIn):
			
			print "CLS: ", cls.__tablename__
			print "Columns: ", queryColumns

			query = session.query(cls.dbid, cls)
			targets = query.filter( or_(qC.ilike(equal+'%') for qC in queryColumns) ).all()
			if targets:
				print "return: ", returnColumns
				return [ (partDBID, geneDBID, {}) for partDBID, geneDBID in zip(
							[ target[0] for target in targets ],
							[ target[1].gene[0].dbid for target in targets ]
						) ]
	
	return {}

def processQuery(sesion, scope, term, displayColumns, nHits):

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
	for scTable in scopeDict:
		scColumns = scopeDict[scTable]
		data = pd.DataFrame(columns=["ID"] + displayColumns)
		newData = queryTable(sesion, scTable, scColumns, term, displayColumns)

		data.append(newData, ignore_index = True, sort = False)

	# 	for hit in newData:
	# 		partDBID = hit[0]
	# 		geneDBID = hit[1]
	# 		cols = { col: '' for col in displayColumns }
	# 		cols.update(hit[2])
	# 		cols["dbid"] = partDBID


	# 		if not geneDBID in genes:
	# 			genes[geneDBID] = []

	# 		genes[geneDBID].append( cols )

	# for geneDBID in genes:
	# 	genes[geneDBID].sort(key = lambda hit: hit["eVal"])

	# sortedKeys = sorted( genes, key = lambda k: genes[k][0]["eVal"] )

	# table = {}
	# table["header"] = ["dbid"] + displayColumns
	# table["data"] = []

	# rowid = 0
	# for geneDBID in sortedKeys:
	# 	rowid += 1
	# 	geneid = rowid
		
	# 	cols = genes[geneDBID][0].copy()
	# 	cols["dbid"] = geneDBID

	# 	row = {"rowid": rowid, "pid": "none", "cols": cols }
	# 	table["data"].append(row)

	# 	for hit in genes[geneDBID][0:5]:
	# 		rowid +=1
	# 		row = {"rowid": rowid, "pid": geneid, "cols": hit }
	# 		table["data"].append(row)

	# for row in table["data"]:
	# 	print row
	return table