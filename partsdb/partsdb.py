from tables import PartMixIn, Base, Sys
from misc import nextIDGenerator

from tools import prepareLibrary, compRev, getUtrCoordinates
from tools import saveSequences

from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr

class PartsDB:

	Session 	= sessionmaker()
	newParts	= []

	def __init__(self, address, Base, clean = False, idGenerator = nextIDGenerator):
		self.engine = create_engine(address, echo = False)
		self.Session.configure(bind=self.engine)
		self.Base = Base

		self.classes = {}

		for cls in Base.__subclasses__():
			if issubclass(cls, PartMixIn):
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

	def populate(self, genomeFileName, transcriptFileName, proteinFileName, mapFileName):
		loci, transcripts = prepareLibrary(genomeFileName, transcriptFileName, proteinFileName, mapFileName)

		session = self.Session()

		n = 0

		for locusName, locus in loci.iteritems():
			locus["PromoterP"] = None
			locus["PromoterN"] = None
			locus["TerminatorP"] = None
			locus["TerminatorN"] = None
			locus["object"] = self.addPart('locus', coordinates = locus["coordinates"] )

		for transcriptName, transcript in transcripts.iteritems():

			# locus = loci[ transcript["locusName"] ]

			transStart = min ( [ex[0] for ex in transcript["exons"]] )
			transStop =  max ( [ex[1] for ex in transcript["exons"]] )

			# print transcriptName, transStart, transStop

			for cdsName, cds in transcript["cdss"].iteritems():
				
				n += 1
				
				if transcript["direction"] == cds["transLoc"][2]:
					direction = '+'
				else:
					direction = '-'

				coordinates = ';'.join(str(cds[0]) + ',' + str(cds[1]) + ',' + direction for cds in cds["geneLoc"])

				cdsStart = min( [ex[0] for ex in  cds["geneLoc"]] )
				cdsStop  = max( [ex[1] for ex in  cds["geneLoc"]] )

				sequence = compRev(transcript["sequence"][cdsStart-1:cdsStop], direction)

				print transcriptName, cdsName, coordinates

				cds["object"] = self.addPart('cds', coordinates= coordinates, seq = sequence)

				if direction == '+':
					if not locus["PromoterP"]:
						sequence = compRev(transcript["sequence"][:transStart-1], direction)

						promoter = self.addPart('promoter', seq = sequence )

						# promoterID = insertRow('promoter', {'seq' : sequence})

						sequence = compRev(transcript["sequence"][transStop-1:], direction)
						
						terminator = self.addPart('terminator', seq = sequence )

						locus["PromoterP"] = promoter
						locus["TerminatorP"] = terminator
					else:
						promoterID = locus["PromoterP"]
						terminatorID = locus["TerminatorP"]
				
				if direction == '-':
					if not locus["PromoterN"]:
						sequence = compRev(transcript["sequence"][transStop-1:], direction)
						
						promoter = self.addPart('promoter', seq = sequence )

						sequence = compRev(transcript["sequence"][:transStart-1], direction)
						
						terminator = self.addPart('terminator', seq = sequence )
						
						locus["PromoterN"] = promoter
						locus["TerminatorN"] = terminator
					else:
						promoterID = locus["PromoterN"]
						terminatorID = locus["TerminatorN"]

				utrLD = {'+':'utr5', '-':'utr3'}
				utrRD = {'+':'utr3', '-':'utr5'}

				sequence = compRev(transcript["sequence"][transStart-1:cdsStart], direction)
				coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)
				
				if direction == '+':
					utr5 = self.addPart('utr5', seq = sequence, coordinates = coordinates )
					# utrL = insertRow(utrLD[direction] , {'seq': sequence, 'coordinates': coordinates})
				else:
					utr3 = self.addPart('utr3', seq = sequence, coordinates = coordinates )

				sequence = compRev(transcript["sequence"][cdsStop-1:transStop], direction)
				coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)

				if direction == '+':
					utr3 = self.addPart('utr3', seq = sequence, coordinates = coordinates )
				else:
					utr5 = self.addPart('utr5', seq = sequence, coordinates = coordinates )


				gene = self.addPart( 'gene', promoter = promoter,  utr5 = utr5, cds = cds["object"], utr3 = utr3, terminator = terminator, locus = locus["object"] )
		print "added {0} cdss".format(n)
		self.commitParts()

	def addPart(self, table, **kwargs):
		newPart = self.classes[table](**kwargs)
		self.newParts.append(newPart)
		return newPart

	def commitParts(self):
		print len(self.newParts)
		if len(self.newParts) != 0:
			session = self.Session()
			session.add_all(self.newParts)
			session.commit()
			prefix = self._getSysVal('prefix')
		
			for clsName, cls in self.classes.iteritems():
				unNamed = session.query(cls).filter(cls.dbid == None).all()
				for obj in unNamed:
					obj.dbid = "{0}.{1}.{2}".format(prefix, clsName, obj.id)
					session.add(obj)
			session.commit()
			session.close()
		self.newParts = []

if __name__ == "__main__":

	genomeFileName = 'data/final.scaffolds.fa'
	transcriptFileName = 'data/trans.fa'
	proteinFileName = 'data/pep.fa'
	mapFileName = 'data/map.gff3'

	# blastFile = prefix + 'output/blastp.info'

	marpodb = PartsDB('postgresql:///testdb', clean = True, Base = Base)
	marpodb.setup(prefix = "mpdb")
	marpodb.addPart('promoter')
	marpodb.commitParts()
	# CDS = marpodb.classes['cds']
	# session = marpodb.Session()
	# session.close()
	# marpodb.populate(genomeFileName, transcriptFileName, proteinFileName, mapFileName)
	# marpodb.annotate('BlastpHit', blastFile)
	# session = marpodb.Session()
	# saveSequences( session.query(UTR5).filter().all(), "cdsSeqs.fasta", translate = False )

