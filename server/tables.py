from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text, Float
from sqlalchemy.orm import relationship
from partsdb.system.Tables import Base, BaseMixIn, PartMixIn, ExonMixIn, AnnotationMixIn
from partsdb.tools.Annotators import BlastAnnotator, PfamAnnotator

class Locus(Base, BaseMixIn):
	coordinates 	= Column( Text )
	genes 			= relationship("Gene", back_populates="locus")

class Promoter(Base,BaseMixIn,PartMixIn):
	gene 			= relationship("Gene", back_populates="promoter")
	geneID 			= Column(Integer, ForeignKey('gene.id'))

class UTR5(Base,BaseMixIn,PartMixIn, ExonMixIn):
	gene 			= relationship("Gene", back_populates="utr5")
	geneID 			= Column(Integer, ForeignKey('gene.id'))

class CDS(Base,BaseMixIn,PartMixIn, ExonMixIn):
	gene 			= relationship("Gene", back_populates="cds")
	geneID 			= Column(Integer, ForeignKey('gene.id'))

class UTR3(Base,BaseMixIn,PartMixIn, ExonMixIn):
	gene 			= relationship("Gene", back_populates="utr3")
	geneID 			= Column(Integer, ForeignKey('gene.id'))

class Terminator(Base,BaseMixIn,PartMixIn):
	gene 			= relationship("Gene", back_populates="terminator")
	geneID 			= Column(Integer, ForeignKey('gene.id'))

class Gene(Base,BaseMixIn):
	name 			= Column( String(100) )
	alias			= Column( String(100) )
	transcriptName  = Column( String(200) )
	# promoterID  	= Column( Integer, ForeignKey('promoter.id') )
	# utr5ID  		= Column( Integer, ForeignKey('utr5.id') )
	# cdsID  			= Column( Integer, ForeignKey('cds.id') )
	# utr3ID  		= Column( Integer, ForeignKey('utr3.id') )
	# terminatorID  	= Column( Integer, ForeignKey('terminator.id') )
	locusID  		= Column( Integer, ForeignKey('locus.id') )
	locusStrand     = Column( Integer )

	promoter 		= relationship(Promoter, 	uselist=False, back_populates="gene")
	utr5 			= relationship(UTR5, 		uselist=False, back_populates="gene")
	cds 			= relationship(CDS, 		uselist=False, back_populates="gene")
	utr3 			= relationship(UTR3, 		uselist=False, back_populates="gene")
	terminator 		= relationship(Terminator, 	uselist=False, back_populates="gene")
	locus   		= relationship(Locus,		uselist=False, back_populates="genes")

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