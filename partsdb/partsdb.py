from tables import Promoter, UTR5, CDS, UTR3, Terminator, Gene, Locus, Sys, Base, BaseMixIn
from misc import nextIDGenerator

from tools import prepareLibrary, compRev, getUtrCoordinates


from sqlalchemy.orm import sessionmaker, mapper
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr
class PartsDB:

	Session 	= sessionmaker()
	metadata 	= MetaData()
	newParts	= []

	def __init__(self, address, prefix='n', classes = { 'promoter' : Promoter, 'utr5' : UTR5,  'cds': CDS, 'utr3' : UTR3, 'terminator' : Terminator, 'gene' : Gene, 'locus' : Locus}, idGenerator = nextIDGenerator):
		self.engine = create_engine(address, echo = False)
		self.Session.configure(bind=self.engine)

		self.classes = {}

		self.idGenerator = idGenerator

		for cls in classes:
			self.classes[cls] = self._partClass( classes[cls], address, self.idGenerator )
			classes[cls].metadata.create_all(self.engine, checkfirst=True)

		Sys.metadata.create_all(self.engine, checkfirst=True)

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

	# def clear(self):
	# 	for cls in self.classes:
	# 		self.classes[cls].metadata.drop_all(self.engine, checkfirst=True)

	def _partClass(self, cls, address, idGenerator):

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

				prefix = session.query(Sys.value).filter(Sys.variable == 'prefix').first()[0]

				self.id = '{0}.{1}.{2}'.format(prefix, cls.__tablename__, nid)

		table = wrapperClass.__table__

		mapper(wrapperClass, table, non_primary = True)

		return(wrapperClass)

	def populate(self, genomeFileName, transcriptFileName, proteinFileName, mapFileName):
		loci, transcripts = prepareLibrary(genomeFileName, transcriptFileName, proteinFileName, mapFileName)

		print loci

		session = self.Session()

		for locusName, locus in loci.iteritems():
			locus["PromoterP"] = None
			locus["PromoterN"] = None
			locus["TerminatorP"] = None
			locus["TerminatorN"] = None
			locus["object"] = self.addPart('locus', coordinates = locus["coordinates"] )

		for transcriptName, transcript in transcripts.iteritems():

			locus = loci[ transcript["locusName"] ]

			transStart = min ( [ex[0] for ex in transcript["exons"]] )
			transStop =  max ( [ex[1] for ex in transcript["exons"]] )


		for cdsName, cds in transcript["cdss"].iteritems():
			if transcript["direction"] == cds["transLoc"][2]:
				direction = '+'
			else:
				direction = '-'

			coordinates = ';'.join(str(cds[0]) + ',' + str(cds[1]) + ',' + direction for cds in cds["geneLoc"])

			cdsStart = min( [ex[0] for ex in  cds["geneLoc"]] )
			cdsStop  = max( [ex[1] for ex in  cds["geneLoc"]] )

			sequence = compRev(locus["seq"][cdsStart-1:cdsStop], direction)

			cds["object"] = self.addPart('cds', coordinates= coordinates, seq = sequence)

			if direction == '+':
				if not locus["PromoterP"]:
					sequence = compRev(locus["seq"][:transStart-1], direction)

					promoter = self.addPart('promoter', seq = sequence )

					# promoterID = insertRow('promoter', {'seq' : sequence})

					sequence = compRev(locus["seq"][transStop-1:], direction)
					
					terminator = self.addPart('terminator', seq = sequence )

					locus["PromoterP"] = promoter
					locus["TerminatorP"] = terminator
				else:
					promoterID = locus["PromoterP"]
					terminatorID = locus["TerminatorP"]
			
			if direction == '-':
				if not locus["PromoterN"]:
					sequence = compRev(locus["seq"][transStop-1:], direction)
					
					promoter = self.addPart('promoter', seq = sequence )

					sequence = compRev(locus["seq"][:transStart-1], direction)
					
					terminator = self.addPart('terminator', seq = sequence )
					
					locus["PromoterN"] = promoter
					locus["TerminatorN"] = terminator
				else:
					promoterID = locus["PromoterN"]
					terminatorID = locus["TerminatorN"]

			utrLD = {'+':'utr5', '-':'utr3'}
			utrRD = {'+':'utr3', '-':'utr5'}

			sequence = compRev(locus["seq"][transStart-1:cdsStart], direction)
			coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)
			
			if direction == '+':
				utr5 = self.addPart('utr5', seq = sequence, coordinates = coordinates )
				# utrL = insertRow(utrLD[direction] , {'seq': sequence, 'coordinates': coordinates})
			else:
				utr3 = self.addPart('utr3', seq = sequence, coordinates = coordinates )

			sequence = compRev(locus["seq"][cdsStop-1:transStop], direction)
			coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)

			if direction == '+':
				utr3 = self.addPart('utr3', seq = sequence, coordinates = coordinates )
			else:
				utr5 = self.addPart('utr5', seq = sequence, coordinates = coordinates )


			gene = self.addPart( 'gene', promoter = promoter,  utr5 = utr5, cds = cds["object"], utr3 = utr3, terminator = terminator, locus = locus["object"] )


	def addPart(self, table, **kwargs):
		session = self.Session()

		newPart = self.classes[table](**kwargs)

		self.newParts.append(newPart)

		session.add(newPart)
		session.commit()
		session.close()
		return newPart

if __name__ == "__main__":

	genomeFileName = '/Users/md/ownCloud/Projects/MarpoDB/marpodb_tools/data/final.scaffolds.fa'
	transcriptFileName = '/Users/md/ownCloud/Projects/MarpoDB/marpodb_tools/output/extra.fa_final.scaffolds.fa/extra.fa.transdecoder.complete.mapped.trans'
	proteinFileName = '/Users/md/ownCloud/Projects/MarpoDB/marpodb_tools/output/extra.fa_final.scaffolds.fa/extra.fa.transdecoder.complete.mapped.pep'
	mapFileName = '/Users/md/ownCloud/Projects/MarpoDB/marpodb_tools/output/extra.fa_final.scaffolds.fa/extra.fa_final.scaffolds.fa.splign.gff3'
	marpodb = PartsDB('postgresql:///testdb')
	marpodb.setup(prefix = "testdb")

	marpodb.populate(genomeFileName, transcriptFileName, proteinFileName, mapFileName)

