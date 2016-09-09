from tables import Promoter, UTR5, CDS, UTR3, Terminator, Gene, Base, BaseMixIn
from misc import nonRemoveIDGenerator

from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr
class PartsDB:

	Session 	= sessionmaker()
	metadata 	= MetaData()

	def __init__(self, address, prefix, classes = { 'promoter' : Promoter, 'utr5' : UTR5,  'cds': CDS, 'utr3' : UTR3, 'terminator' : Terminator, 'gene' : Gene}, idGenerator = nonRemoveIDGenerator):
		self.engine = create_engine(address, echo = True)
		self.Session.configure(bind=self.engine)
		self.prefix = prefix

		self.classes = {}

		self.idGenerator = idGenerator

		for cls in classes:
			self.classes[cls] = self._partClass( classes[cls], address, prefix, self.idGenerator )
			classes[cls].metadata.create_all(self.engine)

		

	def _partClass(self, cls, address, prefix, idGenerator):

		table = cls.__table__

		class wrapperClass(cls):
			@declared_attr
			def __tablename__(clss):
				return cls.__name__.lower()

			def __init__(self, **kwargs):
				super(cls, self).__init__(**kwargs)

				engine = create_engine(address)
				session = sessionmaker(bind=engine)()
				# self.id = idGenerator(cls, session)
				nid = idGenerator(cls, session)


				self.id = '{0}.{1}.{2}'.format(prefix, cls.__tablename__, nid)

		mapper(wrapperClass, table, non_primary = True)

		return(wrapperClass)

	def addPart(self, table, **kwargs):
		session = self.Session()

		newPart = self.classes[table](**kwargs)

		session.add(newPart)
		session.commit()

if __name__ == "__main__":
	marpodb = PartsDB('postgresql:///testdb', 'testdb')
	marpodb.addPart('cds', seq = 'MarpoDBCDS')
	marpodb.addPart('promoter', seq = 'MarpoDBCDS')
	marpodb.addPart('gene', cdsID = 'testdb.cds.0', promoterID = 'testdb.promoter.0')




