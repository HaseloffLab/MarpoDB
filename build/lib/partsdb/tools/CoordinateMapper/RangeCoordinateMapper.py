from CoordinateMapper import CoordinateMapper
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

class RangeCoordinateMapper(CoordinateMapper):
	def __init__(self, selist, length, startOffset, endOffset):
		if selist.strand == -1:
			selist = selist._flip(length)
			self.strand = -1
		else:
			self.strand = 1
		super(RangeCoordinateMapper, self).__init__(selist)
		self.length = length
		self.startOffset = startOffset
		self.endOffset = endOffset
		
	# start - end in GenBank notation
	def rc2g(self, start, end, strand):
		locations = []
		# print self._exons, len(self._exons)
		# print start, end, self.startOffset, self.endOffset, self.strand, strand
		
		start = self.c2g(start - self.startOffset, dialect = 'GenBank').to_genbank()
		end   = self.c2g(end - self.startOffset, dialect = 'GenBank').to_genbank()
		
		# if strand == self.strand:
		# 	start = self.c2g(start - self.startOffset, dialect = 'GenBank').to_genbank()
		# 	end   = self.c2g(end - self.startOffset, dialect = 'GenBank').to_genbank()
		# else:
		# 	start = self.c2g(start - self.endOffset, dialect = 'GenBank').to_genbank()
		# 	end   = self.c2g(end - self.endOffset, dialect = 'GenBank').to_genbank()
		# print start, end


		for exon in sorted(self._exons.parts, key = lambda x: x.start):
			eStart = int(exon.start) + 1
			eEnd   = int(exon.end)

			if start <= eStart:
				if end >= eStart and end <= exon.end:
					locations.append( FeatureLocation(eStart-1, end,strand) )
					break
				elif end > eEnd:
					locations.append( FeatureLocation(eStart-1, exon.end, strand) )
			
			elif start <= eEnd:
				if end <= eEnd:
					locations.append( FeatureLocation(start-1, end, strand) )
					break
				else:
					locations.append( FeatureLocation(start-1, eEnd, strand) )
		
		if strand == 1:
			locations.sort(key = lambda x: x.start)
		if strand == -1:
			locations.sort(key = lambda x: x.start, reverse = True)

		if len(locations) == 1:
			location = locations[0]
		else:
			location = CompoundLocation(locations)

		if self.strand == -1:
			return location._flip(self.length)

		return location

if __name__ == "__main__":
	# {300, 399, -}, {100, 199, -}
	exons = SeqFeature( type = 'mRNA', location =  CompoundLocation( [ FeatureLocation(299, 399,-1), FeatureLocation(99, 199,-1)] ) )
	# { 100, 199 }, {300, 399}
	# exons = SeqFeature( type = 'mRNA', location =  CompoundLocation( [ FeatureLocation(99, 199, 1), FeatureLocation(299, 399, 1)] ) )
	rcm = RangeCoordinateMapper( exons, 399 , 0, 0)
	print rcm.rc2g(50, 180, -1)
