from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float
from sqlalchemy.orm import relationship, backref
from partsdb.system.Tables import Base, BaseMixIn, PartMixIn, ExonMixIn, AnnotationMixIn
from partsdb.tools.Annotators import BlastAnnotator, TableAnnotator

class Locus(Base, BaseMixIn):
	coordinates 	= Column( Text )

class Promoter(Base,BaseMixIn,PartMixIn):
	pass

class UTR5(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class CDS(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class UTR3(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class Terminator(Base,BaseMixIn,PartMixIn):
	pass

class Gene(Base,BaseMixIn):
	name 			= Column( String(100) )
	alias			= Column( String(100) )
	transcriptName  = Column( String(200) )
	promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	utr5ID  		= Column( Integer, ForeignKey('utr5.id') )
	cdsID  			= Column( Integer, ForeignKey('cds.id') )
	utr3ID  		= Column( Integer, ForeignKey('utr3.id') )
	terminatorID  	= Column( Integer, ForeignKey('terminator.id') )
	locusID  		= Column( Integer, ForeignKey('locus.id') )
	locusStrand     = Column( Integer )

	promoter 		= relationship(Promoter,	backref=backref("gene"),	enable_typechecks=False)
	utr5 			= relationship(UTR5,		backref=backref("gene"),	enable_typechecks=False)
	cds 			= relationship(CDS, 		backref=backref("gene"),	enable_typechecks=False)
	utr3 			= relationship(UTR3, 		backref=backref("gene"),	enable_typechecks=False)
	terminator 		= relationship(Terminator, 	backref=backref("gene"),	enable_typechecks=False)
	locus   		= relationship(Locus,		backref=backref("gene"),	enable_typechecks=False)

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
	__annotatorclass__ = TableAnnotator

class DbxRef( Base, BaseMixIn, AnnotationMixIn ):
	__targetclass__	= CDS
	__annotatorclass__ = TableAnnotator
