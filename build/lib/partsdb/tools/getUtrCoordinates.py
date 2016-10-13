def getUtrCoordinates(exons, cdsStart, cdsStop):
		coordinates = []

		for exon in sorted(exons, key = lambda x: x[0]):
			if cdsStart < exon[0]:
				if cdsStop <= exon[1] and cdsStop >= exon[0]:
					coordinates.append( (exon[0], cdsStop) )
					break
				elif cdsStop > exon[1]:
					coordinates.append( exon )
			elif cdsStart <= exon[1]:
				if cdsStop <= exon[1]:
					coordinates.append( (cdsStart, cdsStop) )
					break
				else:
					coordinates.append( (cdsStart, exon[1]) )

		return coordinates