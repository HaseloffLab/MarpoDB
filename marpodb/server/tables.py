from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float
from sqlalchemy.orm import relationship
from partsdb.system.Tables import Base, BaseMixIn, PartMixIn, ExonMixIn, AnnotationMixIn
from partsdb.tools.Annotators import BlastAnnotator, PfamAnnotator

class Locus(Base, BaseMixIn):
	coordinates 	= Column( Text )

class Promoter(Base,BaseMixIn,PartMixIn):
	# locusID 		= Column( Integer, ForeignKey('locus.id') )
	# locus 			= relationship(Locus, 		enable_typechecks=False)

	# locusStrand		= Column( Integer )
	pass

class UTR5(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class CDS(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class UTR3(Base,BaseMixIn,PartMixIn, ExonMixIn):
	pass

class Terminator(Base,BaseMixIn,PartMixIn):
	# locusID 		= Column( Integer, ForeignKey('locus.id') )
	# locus 			= relationship(Locus, 		enable_typechecks=False)

	# locusStrand		= Column( Integer )
	pass

class Gene(Base,BaseMixIn):
	name 			= Column( String(100) )
	alias			= Column( String(100) )
	promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	utr5ID  		= Column( Integer, ForeignKey('utr5.id') )
	cdsID  			= Column( Integer, ForeignKey('cds.id') )
	utr3ID  		= Column( Integer, ForeignKey('utr3.id') )
	terminatorID  	= Column( Integer, ForeignKey('terminator.id') )
	locusID  		= Column( Integer, ForeignKey('locus.id') )
	locusStrand     = Column( Integer )

	promoter 		= relationship(Promoter, 	enable_typechecks=False)
	utr5 			= relationship(UTR5, 		enable_typechecks=False)
	cds 			= relationship(CDS, 		enable_typechecks=False)
	utr3 			= relationship(UTR3, 		enable_typechecks=False)
	terminator 		= relationship(Terminator, 	enable_typechecks=False)
	locus   		= relationship(Locus,		enable_typechecks=False)

class BlastpHit(Base, BaseMixIn, AnnotationMixIn):

	__targetclass__ = CDS
	__annotatorclass__ = BlastAnnotator

	uniID 			= Column( String(100) )
	coverage		= Column( Float )
	qLen			= Column( Integer )
	tLen			= Column( Integer )
	coordinates		= Column( Text )
	eVal 			= Column( Float )
	proteinName 	= Column( Text )
	geneName 		= Column( Text )
	origin 			= Column( Text )

class PfamHit(Base, BaseMixIn, AnnotationMixIn):

	__targetclass__ = CDS
	__annotatorclass__ = PfamAnnotator

	name 			= Column( String(100) )
	acc 			= Column( String(100) )
	eVal 			= Column( Float )
	cVal 			= Column( Float )
	description 	= Column( Text )
	coordinates		= Column( Text )