from . import Exporter
from sqlalchemy.inspection import inspect
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ...system.Tables import PartMixIn, ExonMixIn
from Bio.Alphabet import generic_dna
class GenBankExporter(Exporter):

	def coordinatesToLocation(self, coordinates):
		locationParts = [ FeatureLocation(int(p[0]), int(p[1]), int(p[2]) ) for p in [ s.split(',') for s in coordinates.split(';')] ]
		if len(locationParts) == 1:
			return locationParts[0]
		elif len(locationParts) > 1:
			return CompoundLocation(locationParts)
		else:
			return None

	def export(self, gene, outputFileName=None):
		self.keys = ['promoter' , 'utr5', 'cds', 'utr3', 'terminator']
		parts = [  getattr(gene, key) for key in self.keys ]
		
		strand = gene.locusStrand

		if outputFileName:
			gene = SeqRecord(id = str(gene.dbid).replace('.',''), name = str(gene.name), seq = '' )
		else:
			gene = SeqRecord(id = gene.dbid, name = str(gene.name), seq = '' )

		gene.annotations = {'strand' : strand}

		for partType, part in zip( self.keys, parts):
			l = len(gene)
			if isinstance(part, PartMixIn):
				if isinstance(part, ExonMixIn):
					feature = SeqFeature( type = partType, location = self.coordinatesToLocation(part.coordinates)._shift( l ), id=part.dbid )
				else:
					feature = SeqFeature( type = partType, location = FeatureLocation( l, l + len(part.seq) ), id=part.dbid )

				gene.seq += Seq(part.seq, generic_dna)
				gene.features.append(feature)
		if outputFileName:
			print 'outputFileName IS'
		else:
			print 'outputFileName IS NONE'

		if outputFileName:
			outputFile = open(outputFileName, 'w')
			SeqIO.write(gene, outputFile, "gb")
		else:
			return gene

if __name__ == "__main__":
	exporter = GenBankExporter(1)
