from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Text
from sqlalchemy.orm import relationship

Base = declarative_base()

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
	def annotator(cls):
		return cls.__annotatorclass__(cls)

	@declared_attr
	def targetID(cls):
		return Column( Integer, ForeignKey('{0}.id'.format(cls.__targetclass__.__tablename__) ) )

	@declared_attr
	def target(cls):
		return relationship( cls.__targetclass__, enable_typechecks=False )

class Sys(Base):
	__tablename__ 	= 'sys'
	
	variable		= Column( Text, primary_key = True )
	value			= Column( Text )