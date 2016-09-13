from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float
from sqlalchemy.orm import relationship

from tools import annotateBlastp

class BaseMixIn(object):

	id 		= Column( Integer, 		primary_key = True)
	dbid 	= Column( String(30), 	unique		= True)

	@declared_attr
	def __tablename__(cls):
		return cls.__name__.lower()

class PartMixIn(object):
	seq = Column( Text )

class ExonMixIn(object):
	coordinates = Column( Text )

class AnnotationMixIn(object):

	@declared_attr
	def targetID(cls):
		return Column( Integer, ForeignKey('{0}.id'.format(cls.__name__.lower())) )

	@declared_attr
	def target(cls):
		return relationship( cls.__targetclass__, enable_typechecks=False )

Base = declarative_base()

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
	promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	utr5ID  		= Column( Integer, ForeignKey('utr5.id') )
	cdsID  			= Column( Integer, ForeignKey('cds.id') )
	utr3ID  		= Column( Integer, ForeignKey('utr3.id') )
	terminatorID  	= Column( Integer, ForeignKey('terminator.id') )
	locusID  		= Column( Integer, ForeignKey('locus.id') )

	promoter 		= relationship(Promoter, 	enable_typechecks=False)
	utr5 			= relationship(UTR5, 		enable_typechecks=False)
	cds 			= relationship(CDS, 		enable_typechecks=False)
	utr3 			= relationship(UTR3, 		enable_typechecks=False)
	terminator 		= relationship(Terminator, 	enable_typechecks=False)
	locus   		= relationship(Locus,		enable_typechecks=False)

class BlastpHit(Base, BaseMixIn, AnnotationMixIn):

	__targetclass__ = CDS

	cdsID			= Column( Integer, ForeignKey('cds.id') )
	uniID 			= Column( String(100) )
	coverage		= Column( Float )
	qLen			= Column( Integer )
	tLen			= Column( Integer )
	coordinates		= Column( Text )
	eVal 			= Column( Float )
	proteinName 	= Column( Text )
	geneName 		= Column( Text )
	origin 			= Column( Text )

	@staticmethod
	def annotate(fileName, session):
		annotateBlastp(BlastpHit, fileName, session)

class PfamHit(Base, BaseMixIn, AnnotationMixIn):

	__targetclass__ = CDS

	cdsID			= Column( Integer, ForeignKey('cds.id') )
	name 			= Column( String(100) )
	acc 			= Column( String(100) )
	eVal 			= Column( Float )
	cVal 			= Column( Float )
	description 	= Column( Text )
	coordinates		= Column( Text )

class Sys(Base):
	__tablename__ 	= 'sys'
	
	variable		= Column( Text, primary_key = True )
	value			= Column( Text )