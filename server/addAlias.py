from sqlalchemy.sql.expression import bindparam
from partsdb.partsdb import PartsDB
from tables import Base, Gene

import sys

marpodb = PartsDB('postgresql:///' + sys.argv[1], Base = Base)

data = []

inFile = open( 'data/takblast.outfmt6' )

for line in inFile:
	tabs = line.split()
	dbid = tabs[0]
	phtz = tabs[1]

	alias = phtz.split('.')[0]

	cov = 100*(float(tabs[5])-float(tabs[4])) / float(tabs[2])
	hid = float(tabs[8])

	lq  = tabs[2]
	ls  = tabs[3]
	if cov > 95 and hid > 95 and lq == ls and phtz.split('.')[1] == '1':
		row = {'id' : dbid.split('.')[-1], 'alias' : phtz}
		data.append(row)

print len(data)

session = marpodb.Session()
session.bulk_update_mappings(Gene, data)
session.commit()
session.close()