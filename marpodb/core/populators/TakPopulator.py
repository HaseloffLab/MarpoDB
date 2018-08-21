from ..tables import *

from partsdb.tools.Populators.Populator import Populator

from Bio import SeqIO
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

class TakPopulator(Populator):

	def locationToCoordinates(self, location):

		if location.strand == 1:
			location = location._shift(-location.start)
		elif location.strand == -1:
			location = location._flip(location.end)

		if isinstance(location, CompoundLocation):
			return ";".join( [ "{0},{1},1".format( part.start, part.end ) for part in sorted(location.parts, key = lambda part: part.start) ] )

		elif isinstance(location, FeatureLocation):
			return  "{0},{1},1".format( location.start, location.end )

	def extractSequence(self, scaffoldID, location):
		extractLocation = FeatureLocation( location.start, location.end, location.strand )
		seq = self.scaffoldDict[scaffoldID]

		return str(extractLocation.extract(seq).seq).upper()


	def populate(self, mapFileName, genomeFileName, datasetName, datasetVersion):
		mapFile = open(mapFileName)
		mapRecords = GFF.parse(mapFile)

		dataset = self.db.session.query(Dataset).filter(Dataset.name == datasetName, Dataset.version == datasetVersion).first()
		if not dataset:
			dataset = self.db.addPart("dataset", name = datasetName, version = datasetVersion)

		geneMap = {}

		gene = None

		for record in mapRecords:
			for geneFeature in record.features:
				mRNAFeature = geneFeature.sub_features[0]
			
				if mRNAFeature.qualifiers["ID"][0].split('.')[1] != '1':
					continue

				gene = {"seqID":record.id, "location": mRNAFeature.location, "cds": [], 'utr5': [], 'utr3': []}
				geneMap[geneFeature.qualifiers["Name"][0]] = gene
				
				for feature in mRNAFeature.sub_features:
					if feature.type=="three_prime_UTR":
						gene["utr3"].append(feature.location)

					elif feature.type=="five_prime_UTR":
						gene["utr5"].append(feature.location)

					elif feature.type=="CDS":
						gene["cds"].append(feature.location)
				
				if gene["utr3"]:
					if len(gene["utr3"]) > 1:
						gene["utr3"] = CompoundLocation(gene["utr3"])
					else:
						gene["utr3"] = gene["utr3"][0]

				if gene["utr5"]:
					if len(gene["utr5"]) > 1:
						gene["utr5"] = CompoundLocation(gene["utr5"])
					else:
						gene["utr5"] = gene["utr5"][0]
				
				if len(gene["cds"]) > 1:
					gene["cds"] = CompoundLocation(gene["cds"])
				else:
					gene["cds"] = gene["cds"][0]

		print "\tProcessed {0} genes".format(len(geneMap))

		genomeFile = open(genomeFileName)
		self.scaffoldDict = SeqIO.to_dict( SeqIO.parse(genomeFile, "fasta") )

		n = 0
		m = 0
		for geneName in geneMap:
			utr3 = None
			utr5 = None
			geneRec = geneMap[geneName]
			seqID = geneRec["seqID"]
			mRNAStart = geneRec["location"].start
			mRNAEnd = geneRec["location"].end
			mRNAStrand = geneRec["location"].strand

			pepSeq = geneRec["cds"].extract( self.scaffoldDict[geneRec["seqID"]] ).seq.translate()
			if not ( pepSeq.startswith('M') and pepSeq.find('*')==len(pepSeq)-1 ):
				n += 1

			leftBorder = max( mRNAStart-3000, 0 )
			rightBorder = min( mRNAEnd+3000, len(self.scaffoldDict[seqID]) )

			if mRNAStrand == 1:
				promoterLocation = FeatureLocation( leftBorder, max( mRNAStart-1, 0), 1 )
				terminatorLocation = FeatureLocation( mRNAEnd+1, rightBorder, 1 )
			else:
				promoterLocation = FeatureLocation( mRNAEnd+1, rightBorder, -1 )
				terminatorLocation = FeatureLocation( leftBorder, max( mRNAStart-1, 0), -1 )

			promoter = self.db.addPart("promoter", seq = self.extractSequence(seqID, promoterLocation) )
			terminator = self.db.addPart("terminator", seq = self.extractSequence(seqID, terminatorLocation) )

			if geneRec["utr5"]:
				utr5 = self.db.addPart("utr5", seq = self.extractSequence(seqID, geneRec["utr5"]), coordinates = self.locationToCoordinates( geneRec["utr5"] ) )

			if geneRec["utr3"]:
				utr3 = self.db.addPart("utr3", seq = self.extractSequence(seqID, geneRec["utr3"]), coordinates = self.locationToCoordinates( geneRec["utr3"] ) )

			extendedGeneLocation = FeatureLocation( leftBorder, rightBorder, geneRec["location"].strand )

			cds = self.db.addPart("cds", seq = self.extractSequence( seqID, geneRec["cds"] ), coordinates = self.locationToCoordinates(geneRec["cds"]) )

			dbref = self.db.addPart("dbxref", description = geneName, refID = geneName, origin = "Phytozome MP-Tak 3.1", target = cds)

			locus = self.db.addPart("locus", coordinates = "{0} {1}:{2}:{3}".format(seqID, mRNAStart, mRNAEnd, mRNAStrand))

			gene = self.db.addPart("gene", promoter = promoter, terminator = terminator, cds = cds, utr5 = utr5, utr3=utr3, locus = locus, dataset = dataset)
			m += 1
		print("\tAdded {0} genes, {1} had bad CDS".format(m,n))	
		self.db.commit()