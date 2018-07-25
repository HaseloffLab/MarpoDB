from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float, UniqueConstraint
from sqlalchemy.orm import relationship, backref
from partsdb.system.Tables import Base, BaseMixIn, ExonMixIn, AnnotationMixIn
from partsdb.tools.Annotators import TableAnnotator
from sqlalchemy.ext.hybrid import hybrid_property

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from .annotators import BlastAnnotator, InterproAnnotator

class PartMixIn(object):
	seq = Column( Text )

	@hybrid_property
	def record(self):
		return SeqRecord( seq = Seq(self.seq), id = self.dbid, description = "" )


def coordinatesToLocation(coordinates):
		locationParts = [ FeatureLocation(int(p[0]), int(p[1]), int(p[2]) ) for p in [ s.split(',') for s in coordinates.split(';')] ]
		if len(locationParts) == 1:
			return locationParts[0]
		elif len(locationParts) > 1:
			return CompoundLocation(locationParts)
		else:
			return None

class Dataset(Base, BaseMixIn):
	name	= Column( String(100) )
	version = Column( String(100) )
	__table_args__ = (UniqueConstraint('name', 'version', name='dataset_fullname'), )

class Locus(Base, BaseMixIn):
	coordinates 	= Column( Text )

class Promoter(Base,BaseMixIn,PartMixIn):
	pass

class UTR5(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class CDS(Base,BaseMixIn,PartMixIn, ExonMixIn):
	
	@hybrid_property
	def pepRecord(self):
		record = self.record
		
		location = coordinatesToLocation(self.coordinates)
		record.seq = location.extract(record.seq).translate()

		return record

class UTR3(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class Terminator(Base,BaseMixIn,PartMixIn):
	pass

class Gene(Base,BaseMixIn):
	name 			= Column( String(100) )
	transcriptName  = Column( String(200) )
	promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	utr5ID  		= Column( Integer, ForeignKey('utr5.id') )
	cdsID  			= Column( Integer, ForeignKey('cds.id') )
	utr3ID  		= Column( Integer, ForeignKey('utr3.id') )
	terminatorID  	= Column( Integer, ForeignKey('terminator.id') )
	locusID  		= Column( Integer, ForeignKey('locus.id') )
	datasetID		= Column( Integer, ForeignKey('dataset.id') )
	locusStrand     = Column( Integer )

	promoter 		= relationship(Promoter,	backref=backref("gene"),	enable_typechecks=False)
	utr5 			= relationship(UTR5,		backref=backref("gene"),	enable_typechecks=False)
	cds 			= relationship(CDS, 		backref=backref("gene"),	enable_typechecks=False)
	utr3 			= relationship(UTR3, 		backref=backref("gene"),	enable_typechecks=False)
	terminator 		= relationship(Terminator, 	backref=backref("gene"),	enable_typechecks=False)
	locus   		= relationship(Locus,		backref=backref("gene"),	enable_typechecks=False)
	dataset			= relationship(Dataset,		backref=backref("gene"),	enable_typechecks=False)

	# @hybrid_property
	# def record(self):
	# 	keys = ["promoter", "utr5", "cds", "utr3", "terminator"]
	# 	parts = [  getattr(self, key) for key in keys ]

	# 	print "### DEBUG"
	# 	print self.dbid
	# 	print self
	# 	print "### DEBUG"

	# 	record = SeqRecord(id = self.dbid, name = str(self.dbid), seq = '' )

	# 	for partType, part in zip( keys, parts):
	# 		l = len(self.seq)
	# 		if isinstance(part, PartMixIn):
	# 			if isinstance(part, ExonMixIn):
	# 				feature = SeqFeature( type = partType, location = coordinatesToLocation(part.coordinates)._shift( l ), id=part.dbid )
	# 			else:
	# 				feature = SeqFeature( type = partType, location = FeatureLocation( l, l + len(part.seq) ), id=part.dbid )

	# 			record.seq += Seq(part.seq, generic_dna)
	# 			record.features.append(feature)

	# 	return record

	@hybrid_property
	def seq(self):
		seq = ''

		if self.promoter:
			seq = seq + self.promoter.seq
		
		if self.utr5:
			seq = seq + self.utr5.seq
		
		seq = seq + self.cds.seq
		
		if self.utr3:
			seq = seq + self.utr3.seq

		if self.terminator:
			seq = seq + self.terminator.seq

		return seq

class BlastpHit( Base, BaseMixIn, AnnotationMixIn):
	__targetclass__	= CDS
	__annotatorclass__ = BlastAnnotator
	coverage		= Column( Float )
	qLen			= Column( Integer )
	tLen			= Column( Integer )
	tStart			= Column( Integer )
	tEnd 			= Column( Integer )
	identity		= Column( Float )

class InterProHit( Base, BaseMixIn, AnnotationMixIn):
	__targetclass__	= CDS
	__annotatorclass__ = InterproAnnotator

class DbxRef( Base, BaseMixIn, AnnotationMixIn ):
	__targetclass__	= CDS
	__annotatorclass__ = TableAnnotator
