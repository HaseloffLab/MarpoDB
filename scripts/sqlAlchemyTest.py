from sqlalchemy import Table, MetaData, Column, Integer, String, ForeignKey
from sqlalchemy.orm import mapper

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy import Sequence

Base = declarative_base()


class CDS(Base):
	__tablename__ = 'cds'

	id = Column( String(30), primary_key = True )
	name = Column( String(30) )

class Gene(Base):
	__tablename__ = 'gene'
	id 		= Column( String(30), primary_key = True )
	name 	= Column( String(30) )
	cds_id  = Column( String(30), ForeignKey('cds.id') )

class PartsDB:

	Session 	= sessionmaker()
	metadata 	= MetaData()

	def __init__(self, address, prefix, classes = {'cds': CDS, 'gene' : Gene}):
		self.engine = create_engine(address, echo = True)
		self.Session.configure(bind=self.engine)
		self.prefix = prefix

		self.classes = {}

		for cls in classes:
			self.classes[cls] = self._partClass( classes[cls], address, prefix )
			classes[cls].metadata.create_all(self.engine)

	def _partClass(self, cls, address, prefix):

		table = cls.__table__

		class wrapperClass(cls):
			def __init__(self, **kwargs):
				super(cls, self).__init__(**kwargs)

				engine = create_engine(address)
				session = sessionmaker(bind=engine)()

				nid = session.query(cls).count()


				self.id = '{0}.{1}.{2}'.format(prefix, cls.__tablename__, nid)

		mapper(wrapperClass, table, non_primary = True)

		return(wrapperClass)

	def addPart(self, table, **kwargs):
		session = self.Session()

		newPart = self.classes[table](**kwargs)

		session.add(newPart)
		session.commit()


marpodb = PartsDB('postgresql:///testdb', 'testdb')
marpodb.addPart('cds', name = 'MarpoDBCDS')
marpodb.addPart('gene', name = 'FirstGene', cds_id = 'testdb.cds.6')

# engine = create_engine('postgresql:///testdb', echo = True)
# Session = sessionmaker(bind=engine)

# session = Session()
# metadata = MetaData()

# cds = partClass(CDS, 'postgresql:///testdb', 'testdb')(name = "WrappedAgain")

# session.add(cds)
# session.commit()



# class PartsDB:

# 	Session = sessionmaker()
# 	metadata 	= MetaData()
	

# 	class PartsBase(object):

# 		@declared_attr
# 		def __tablename__(cls):
# 			return cls.__name__.lower()

		

# 	Base = declarative_base(cls = PartsBase)

# 	class CDS(Base):
# 		id = Column( Integer, primary_key=True )
# 		name = Column( String(30) )
		
	
# 	def __init__(self, address):
# 		self.engine 		= create_engine(address, echo = True)
# 		self.Session.configure(bind=self.engine)
# 		self.session = self.Session()

# 	def createTables(self):
# 		self.Base.metadata.create_all(self.engine)

# 	def addCDS(self, name):
# 		newCDS = self.CDS(name = name)
# 		self.session.add(newCDS)
# 		self.session.commit()
 

# marpodb = PartsDB('postgresql:///testdb')
# # marpodb.createTables()
# marpodb.addCDS(name="newCDS")

# engine = create_engine('postgresql:///marpodb', echo =True)

# Session = sessionmaker(bind=engine)

# metadata = MetaData()

# cds = Table('cds', metadata,
# 			Column('id', String(30), primary_key=True),
# 			Column('name', String(30))
# 	)

# class CDS(object):
# 	def __init__(self, name):
# 		session = Session()
# 		nid = session.query(CDS).count()
# 		self.name = name
# 		self.id = 'cds.{0}'.format(nid + 1)

# mapper(CDS,cds)

# metadata.create_all(engine)

# session = Session()

# cds1 = CDS(name="SECOND")
# session.add(cds1)
# session.commit()
