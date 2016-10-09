from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from pprint import pprint
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

from CoordinateMapper import RangeCoordinateMapper

def prepareLibrary(transMapFileName, transFileName, proteinFileName, genomeFileName):
	
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
			
			if pepStrand == '-':
				pepStrand = -1
			else:
				pepStrand = 1

			if transcriptName in transcripts:
				transcript = transcripts[transcriptName]
				exons = transcript.features[0]

				gene = transcript

				# Annotating CDS
				rcm = RangeCoordinateMapper(exons, len(transcript), transcript.annotations["startOffset"], transcript.annotations["endOffset"] )
				
				location = rcm.rc2g(pepStart, pepEnd, pepStrand)	
				cdsFeature = SeqFeature(type = 'cds', location = location )

				gene.features.append(cdsFeature)

				gene.annotations["strand"] = pepStrand

				# Annotating UTRs
				location = rcm.rc2g(transcript.annotations["startOffset"] + 1, pepStart - 1, pepStrand)
				utrType  = 'utr5' if pepStrand == 1 else 'utr3'
				utrRFeature = SeqFeature(type = utrType, location = location)
				# utrRFeature = SeqFeature(type = 'utrR', location = location)

				location = rcm.rc2g(pepEnd + 1, len(exons) + transcript.annotations["startOffset"], pepStrand)
				utrType  = 'utr5' if pepStrand == -1 else 'utr3'
				utrLFeature = SeqFeature(type = utrType, location = location)
				# utrLFeature = SeqFeature(type = 'utrL', location = location)

				# if cdsFeature.location.strand == exons.location.strand:
				# 	utrRFeature.type = 'utr5'
				# 	utrLFeature.type = 'utr3'
				# else:
				# 	utrRFeature.type = 'utr3'
				# 	utrLFeature.type = 'utr5'

				gene.features.append(utrRFeature)
				gene.features.append(utrLFeature)
				
				# Annotating Promoter / Terminator
				location  = FeatureLocation( 0, exons.location.start, cdsFeature.location.strand ) 
				partType  = 'promoter' if cdsFeature.location.strand == 1 else 'terminator'
				rightPart = SeqFeature(type = partType, location = location)

				location  = FeatureLocation( exons.location.end, len(gene), cdsFeature.location.strand )
				partType  = 'promoter' if cdsFeature.location.strand == -1 else 'terminator'
				leftPart  = SeqFeature(type = partType, location = location)

				gene.features.append(rightPart)
				gene.features.append(leftPart)


				genes[pepName] = gene

	proteinFile.close()

	return genes

def locationToCoordinates(location):
	if isinstance(location, CompoundLocation):
		return ";".join( [ "{0},{1},{2}".format( part.start, part.end, part.strand ) for part in location.parts ] )

	elif isinstance(location, FeatureLocation):
		return  "{0},{1},{2}".format( location.start, location.end, location.strand )


def populate(db, genes):

	session 	= db.Session()
	Locus 		= db.classes['locus']
	Promoter	= db.classes['promoter']
	Terminator	= db.classes['terminator']

	for geneName, gene in genes.iteritems():
		start = min( [pos for pos in gene.features[0].location] ) + 1 # +1 GenBank notation
		stop  = min( [pos for pos in gene.features[0].location] )

		cutStart = max( 1, start - 3000 )
		cutStop  = max( len(gene), stop + 3000 )

		locusCoordinates = "{0}:{1}-{2}".format( gene.annotations["Locus"], cutStart, cutStop )

		gene = gene[cutStart - 1: cutStop]

		locus =  session.query( Locus ).filter( Locus.coordinates == locusCoordinates ).first()
		if not locus:
			locus = db.addPart(session, 'locus', coordinates = locusCoordinates)

		mRNAStart 	= min(gene.features[0].location)
		mRNAStop 	= max(gene.features[0].location)

		parts = {}

		for feature in gene.features:
			if feature.location.strand == 1:
				location = feature.location.shift( -feature.location.start )
				seq = gene.seq[ feature.start : feature.end ]
			else:
				location = feature.location._flip()
				location.shift( -location.start )
				seq = gene.seq[ feature.start : feature.end ].reverse_comlement()
			parts[ feature.type ] = { "seq" : seq, "coordinates" : locationToCoordinates(location) }


		strand = gene.annotations["strand"]

		promoter = session.query(Promoter).filter( Promoter.locusID == locus.id, Promoter.locusStrand == strand ).first()
		if not promoter:
			promoter = db.addPart(session, 'promoter', seq = features['promoter'].extract(gene.seq), locus = locus, locusStrand = strand )

		terminator = session.query(Terminator).filter( Terminator.locusID == locus.id, Terminator.locusStrand == strand ).first()
		if not terminator:
			terminator = db.addPart(session, 'terminator', seq = features['terminator'].extract(gene.seq), locus = locus, locusStrand = strand )


		cds  = db.addPart(session, 'cds', seq = parts['cds']['seq'], coordinates = parts['cds']['coordinates'] )
		utr5 = db.addPart(session, 'utr5', seq = parts['utr5']['seq'], coordinates = parts['utr5']['coordinates'] )
		utr3 = db.addPart(session, 'utr3', seq = parts['utr3']['seq'], coordinates = parts['utr3']['coordinates'] )

		db.addPart(session, 'gene', cds = cds, utr3 = utr3, utr5 = utr5, promoter = promoter, terminator = terminator, locus = locus)
		
if __name__ == '__main__':
	genes = prepareLibrary('/Users/md/marpodb/partsdb/data/map.gff3', '/Users/md/marpodb/partsdb/data/trans.fa', '/Users/md/marpodb/partsdb/data/pep.fa',  '/Users/md/marpodb/partsdb/data/genome.fa')

	for geneName, gene in genes.iteritems():
		print geneName, gene.annotations["Locus"]
		print 'SEQ', gene.seq
		for feature in gene.features:
			print feature
			print feature.extract(gene.seq)

	# populate(1, genes)



