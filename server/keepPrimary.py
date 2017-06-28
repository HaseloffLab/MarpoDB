import sys
from partsdb.partsdb import PartsDB
from tables import *
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation


def coordinatesToLocation(coordinates):
	locationParts = [ FeatureLocation(int(p[0]), int(p[1]), int(p[2]) ) for p in [ s.split(',') for s in coordinates.split(';')] ]
	if len(locationParts) == 1:
		return locationParts[0]
	elif len(locationParts) > 1:
		return CompoundLocation(locationParts)
	else:
		return None


marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

session = marpodb.Session()
n = 0
for locus in session.query(Locus).all():
	genes = session.query(Gene).filter(Gene.locusID == locus.id).all()
	if not genes:
		session.delete(locus)
	else:
		if len(genes) > 1:
			primaryGene = max( genes, key = lambda gene: len( coordinatesToLocation(gene.cds.coordinates) ) )
			genes.remove(primaryGene)
			for gene in genes:
				session.delete(gene.cds)
				session.delete(gene.utr5)
				session.delete(gene.utr3)
				session.delete(gene)
	session.commit()