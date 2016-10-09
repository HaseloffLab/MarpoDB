from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr

from system.Tables import BaseMixIn, Base, Sys
from system.IDGenerator import nextIDGenerator


class PartsDB:

	Session 	= sessionmaker()
	newParts	= []

	def __init__(self, address, Base, clean = False, idGenerator = nextIDGenerator):
		self.engine = create_engine(address, echo = False)
		self.Session.configure(bind=self.engine)
		self.Base = Base

		self.session = self.Session()

		self.classes = {}

		for cls in Base.__subclasses__():
			if issubclass(cls, BaseMixIn):
				self.classes[cls.__tablename__] = cls

		self.idGenerator = idGenerator

		if clean:
			self.Base.metadata.drop_all(self.engine)

		Base.metadata.create_all(self.engine, checkfirst=True)

	def setup(self, **kwargs):
		session = self.Session()
		
		for key, value in kwargs.iteritems():
			val = session.query(Sys).filter(Sys.variable == key).first()
			if val:
				val.value = value
			else:
				val = Sys( variable = key, value = value )
			session.add(val)
		session.commit()

	def _getSysVal(self, key):
		session = self.Session()
		val = session.query(Sys.value).filter(Sys.variable == key).first()
		return val[0]

	def annotate(self, tableName, fileName):
		session = self.Session()
		self.annotationTables[tableName].annotate(fileName, session)
		session.commit()
		session.close()

	def addPart(self, table, **kwargs):
		newPart = self.classes[table](**kwargs)
		self.session.add(newPart)
		return newPart

	def commit(self):
		self.session.commit()
		prefix = self._getSysVal('prefix')
		
		for clsName, cls in self.classes.iteritems():
			unNamed = self.session.query(cls).filter(cls.dbid == None).all()
			for obj in unNamed:
				obj.dbid = "{0}.{1}.{2}".format(prefix, clsName, obj.id)
				self.session.add(obj)
		self.session.commit()

if __name__ == "__main__":

	from tables import *
	from tools.Populators import PlantPopulator
	from tools.Exporters import GenBankExporter

	genomeFileName = 'data/final.scaffolds.fa'
	transcriptFileName = 'data/trans.fa'
	proteinFileName = 'data/pep.fa'
	mapFileName = 'data/map.gff3'

	# blastFile = prefix + 'output/blastp.info'
	marpodb = PartsDB('postgresql:///testdb', clean = True, Base = Base)
	marpodb.setup(prefix = "mpdb")

	ppl = PlantPopulator(marpodb)
	ppl.populate('/Users/md/marpodb/partsdb/data/map.gff3', '/Users/md/marpodb/partsdb/data/trans.fa', '/Users/md/marpodb/partsdb/data/pep.fa',  '/Users/md/marpodb/partsdb/data/genome.fa')

	session = marpodb.Session()
	genes = session.query(Gene).all()

	exporter = GenBankExporter(marpodb)

	for gene in genes:
		exporter.export(gene, gene.dbid)

	# CDS = marpodb.classes['cds']
	# session = marpodb.Session()
	# session.close()
	# marpodb.populate(genomeFileName, transcriptFileName, proteinFileName, mapFileName)
	# marpodb.annotate('BlastpHit', blastFile)
	# session = marpodb.Session()
	# saveSequences( session.query(UTR5).filter().all(), "cdsSeqs.fasta", translate = False )

