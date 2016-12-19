from . import Populator

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from pprint import pprint
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

from ..CoordinateMapper import RangeCoordinateMapper

class PlantPopulator(Populator):

	def _locationToCoordinates(self, location):
		if isinstance(location, CompoundLocation):
			return ";".join( [ "{0},{1},1".format( part.start, part.end ) for part in location.parts ] )

		elif isinstance(location, FeatureLocation):
			return  "{0},{1},1".format( location.start, location.end )

	def populate(self, *args, **kwargs):
		if isinstance(args[1],list):
			self.populateFromList(args[0])
		elif isinstance(args[1], str):
			# print args
			self.populateFromFile(*args)

	def populateFromFile(self, transMapFileName, transFileName, proteinFileName, genomeFileName):
		transMapFile 	= open(transMapFileName)
		transFile 		= open(transFileName)
		proteinFile   	= open(proteinFileName)
		genomeFile 		= open(genomeFileName)
		
		scaffoldDict = SeqIO.to_dict( SeqIO.parse(genomeFile, "fasta") )
		genomeFile.close()

		transMapRecords = GFF.parse(transMapFile, base_dict = scaffoldDict)

		transcripts = {}

		for record in transMapRecords:
			for feature in record.features:
				target = feature.qualifiers["Target"][0]
				transcriptName  = target.split()[0]

				feature.qualifiers 	= { }
				feature.type 		= "mRNA"
				feature.id 			= transcriptName 
				if not transcriptName in transcripts:
					transcript = SeqRecord( record.seq, annotations = {"Locus" : record.id} )
					transcript.features = [ feature ]
					transcripts[transcriptName] = transcript
					transcripts[transcriptName].annotations["minTarget"] = int(target.split()[1])
					transcripts[transcriptName].annotations["maxTarget"] = int(target.split()[2])
				else:
					transcripts[transcriptName].features[0].location += feature.location
					transcripts[transcriptName].annotations["minTarget"] = min([ transcripts[transcriptName].annotations["minTarget"], int(target.split()[1]) ])
					transcripts[transcriptName].annotations["maxTarget"] = max([ transcripts[transcriptName].annotations["maxTarget"], int(target.split()[2]) ])

		transMapFile.close()

		for line in transFile:
			if line.startswith('>'):
				transcriptName = line.split()[0][1:]
				if transcriptName in transcripts:
					length   = int(line.split()[1].split('=')[1])
					transcripts[transcriptName].annotations["startOffset"] = transcripts[transcriptName].annotations["minTarget"] - 1
					transcripts[transcriptName].annotations["endOffset"]   = length - transcripts[transcriptName].annotations["maxTarget"]

		transFile.close()

		genes = {}

		for line in proteinFile:
			if line.startswith('>'):
				substring = line.split()[-1]
				transcriptName = substring.split(':')[0]
				pepName 	   = line.split()[0][1:]
				pepStart = int(substring.split(':')[1].split('-')[0])
				pepEnd   = int(substring.split(':')[1].split('-')[1].split('(')[0] )
				pepStrand = substring.split('(')[1][0]
				# print 'pepStrand', pepStrand
				
				# print pepName

				if pepStrand == '-':
					pepStrand = -1
				else:
					pepStrand = 1

				if transcriptName in transcripts:
					transcript = transcripts[transcriptName]

					gene = transcript

					#Cutting the gene
					start = gene.features[0].location.start
					end   = gene.features[0].location.end

					cutStart = max( 0, start - 3000 )
					cutEnd  = min( len(gene), end + 3000 )

					gene = gene[cutStart:cutEnd]
					gene.annotations = transcript.annotations
					gene.annotations["LocusCoordinates"] = "{0}:{1}-{2}".format( gene.annotations["Locus"], cutStart, cutEnd )
		
					# Annotating CDS
					exons = gene.features[0]
					# print transcriptName
					
					rcm = RangeCoordinateMapper(exons, len(gene), gene.annotations["startOffset"], gene.annotations["endOffset"] )
					
					try:
						location = rcm.rc2g(pepStart, pepEnd, pepStrand)	
					except:
						continue
					cdsFeature = SeqFeature(type = 'cds', location = location )

					print "Pep: ", pepStart, pepEnd, pepStrand
					print "Location: ", location

					gene.features.append(cdsFeature)


					# Annotating UTRs - OLD
					# if gene.annotations["startOffset"] + 1 <= pepStart - 1:
					# 	location = rcm.rc2g(gene.annotations["startOffset"] + 1, pepStart - 1, cdsFeature.location.strand)
					# 	utrType  = 'utr5' if cdsFeature.location.strand == 1 else 'utr3'
					# 	utrRFeature = SeqFeature(type = utrType, location = location)
					# else:
					# 	utrRFeature = None

					# if pepEnd + 1 <= len(exons) + gene.annotations["startOffset"]:
					# 	location = rcm.rc2g(pepEnd + 1, len(exons) + gene.annotations["startOffset"], cdsFeature.location.strand)
					# 	utrType  = 'utr5' if cdsFeature.location.strand == -1 else 'utr3'
					# 	utrLFeature = SeqFeature(type = utrType, location = location)
					# else:
					# 	utrL = None

					# utrLFeature = SeqFeature(type = 'utrL', location = location)

					# if cdsFeature.location.strand == exons.location.strand:
					# 	utrRFeature.type = 'utr5'
					# 	utrLFeature.type = 'utr3'
					# else:
					# 	utrRFeature.type = 'utr3'
					# 	utrLFeature.type = 'utr5'

					# gene.features.append(utrRFeature)
					# gene.features.append(utrLFeature)

					# Annotating UTRs - NEW

					if pepStart > 1:
						location = rcm.rc2g( 1, pepStart-1, pepStrand )
						utrType = 'utr5' if pepStrand == 1 else 'utr3'
						utrLFeature = SeqFeature( type = utrType, location = location )
					else:
						utrLFeature = None


					if pepEnd < len(exons):
						location = rcm.rc2g( pepEnd+1, len(exons), pepStrand )
						utrType = 'utr3' if pepStrand == 1 else 'utr5'
						utrRFeature = SeqFeature( type = utrType, location = location )
					else:
						utrRFeature = None


					gene.features.append(utrRFeature)
					gene.features.append(utrLFeature)
					
					# Annotating Promoter / Terminator
					location  = FeatureLocation( 0, exons.location.start, cdsFeature.location.strand ) 
					partType  = 'promoter' if cdsFeature.location.strand == 1 else 'terminator'
					rightPart = SeqFeature(type = partType, location = location)

					location  = FeatureLocation( exons.location.end, len(gene), cdsFeature.location.strand )
					partType  = 'terminator' if cdsFeature.location.strand == 1 else 'promoter'
					leftPart  = SeqFeature(type = partType, location = location)

					gene.features.append(rightPart)
					gene.features.append(leftPart)


					genes[pepName] = gene

		proteinFile.close()

		self.populateFromList(genes)

	def populateFromList(self, genes):
		assert 'locus' 		in self.db.classes
		assert 'promoter' 	in self.db.classes
		assert 'terminator' in self.db.classes
		
		session 	= self.db.session
		Locus 		= self.db.classes['locus']
		Promoter	= self.db.classes['promoter']
		Terminator	= self.db.classes['terminator']
		Gene 		= self.db.classes['gene']

		for geneName, gene in genes.iteritems():
			# start = min( [pos for pos in gene.features[0].location] ) + 1 # +1 GenBank notation
			# stop  = min( [pos for pos in gene.features[0].location] )

			# cutStart = max( 1, start - 3000 )
			# cutStop  = max( len(gene), stop + 3000 )

			# locusCoordinates = "{0}:{1}-{2}".format( gene.annotations["Locus"], cutStart, cutStop )

			locusCoordinates = gene.annotations["LocusCoordinates"]

			locus =  session.query( Locus ).filter( Locus.coordinates == locusCoordinates ).first()
			if not locus:
				locus = self.db.addPart('locus', coordinates = locusCoordinates)

			mRNAStart 	= min(gene.features[0].location)
			mRNAStop 	= max(gene.features[0].location)

			parts = {}

			strand = 0

			for feature in gene.features:
				if feature:

					if feature.type == 'cds':
						strand = feature.location.strand

					if feature.location.strand == 1:
						location = feature.location._shift( -feature.location.start )
						seq = str(gene.seq[ feature.location.start : feature.location.end ])
					else:
						location = feature.location._flip( len(gene) )
						location = location._shift( -location.start )
						seq = str(gene.seq[ feature.location.start : feature.location.end ].reverse_complement())
					parts[ feature.type ] = { "seq" : seq, "coordinates" : self._locationToCoordinates(location) }
				# print feature.type, location, seq
			# print parts.keys()
			
			promoter = session.query(Promoter).filter( Promoter.id == Gene.promoterID ).\
								filter(Gene.locusID == locus.id, Gene.locusStrand == strand).first()
			if not promoter:
				promoter = self.db.addPart('promoter', seq = parts['promoter']['seq'])
			
			# promoter = session.query(Promoter).filter( Promoter.locusID == locus.id, Promoter.locusStrand == strand ).first()
			# if not promoter:
			# 	promoter = self.db.addPart('promoter', seq = parts['promoter']['seq'], locus = locus, locusStrand = strand )

			terminator = session.query(Terminator).filter( Terminator.id == Gene.terminatorID ).\
								filter(Gene.locusID == locus.id, Gene.locusStrand == strand).first()
			if not terminator:
				terminator = self.db.addPart('terminator', seq = parts['terminator']['seq'])
			
			# terminator = session.query(Terminator).filter( Terminator.locusID == locus.id, Terminator.locusStrand == strand ).first()
			# if not terminator:
			# 	terminator = self.db.addPart('terminator', seq = parts['terminator']['seq'], locus = locus, locusStrand = strand )


			cds  = self.db.addPart('cds', seq = parts['cds']['seq'], coordinates = parts['cds']['coordinates'] )
			
			newGene = self.db.addPart('gene', cds = cds, promoter = promoter, terminator = terminator, locus = locus, locusStrand = strand)

			if 'utr5' in parts:
				utr5 = self.db.addPart('utr5', seq = parts['utr5']['seq'], coordinates = parts['utr5']['coordinates'] )
				newGene.utr5 = utr5
			if 'utr3' in parts:
				utr3 = self.db.addPart('utr3', seq = parts['utr3']['seq'], coordinates = parts['utr3']['coordinates'] )
				newGene.utr3 = utr3
			print newGene.id

		self.db.commit()	