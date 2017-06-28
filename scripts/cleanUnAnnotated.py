import sys
from partsdb.partsdb import PartsDB
from tables import *

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

session = marpodb.Session()

n=0
i=0

for gene in session.query(Gene).all():
	n = session.query(InterProHit).filter(InterProHit.targetID == gene.cds.id).count()
	
	if n == 0:
		session.query(BlastpHit).filter(BlastpHit.targetID == gene.cds.id).delete()
		session.delete(gene)
		n = n+1

	if i % 10 == 0:
		print "Processed {0} genes"

	i = i+1

session.commit()
print "Deleted {0} genes".format(n)