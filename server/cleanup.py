import sys
from partsdb.partsdb import PartsDB
from tables import *

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

session = marpodb.Session()

CDSs = session.query(CDS).all()

for cds in CDSs:
	nBlast =  session.query(BlastpHit).filter(BlastpHit.targetID == cds.id).count()
	nPfam  =  session.query(PfamHit).filter(PfamHit.targetID == cds.id).count()

	if nBlast + nPfam == 0:
		gene = session.query(Gene).filter(Gene.cdsID == cds.id).first()
		session.delete( cds )
		session.delete( gene )

promoters = session.query(Promoter).all()

for promoter in promoters:
	if session.query(Gene).filter(Gene.promoterID == promoter.id).count() == 0:
		session.delete(promoter)

terminators = session.query(Terminator).all()

for terminator in terminators:
	if session.query(Gene).filter(Gene.terminatorID == terminator.id).count() == 0:
		session.delete(terminator)

session.commit()
session.close()


