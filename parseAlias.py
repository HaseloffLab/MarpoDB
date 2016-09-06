import sys
import psycopg2

inFile = open(sys.argv[1])

conn = psycopg2.connect("dbname=marpodb")
cur = conn.cursor()

for line in inFile:
	tabs = line.split()
	camName = tabs[0]
	alias = tabs[1]

	cur.execute("UPDATE gene SET alias = %s WHERE name = %s RETURNING id", (alias, camName))
	res = cur.fetchone()
	if res == None:
		print "{0} not found".format(camName)

conn.commit()
cur.close()
conn.close()
