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

	genes = {}

	for record in transMapRecords:
		for feature in record.features:
			target = feature.qualifiers["Target"][0]
			geneName  = target.split()[0]

			feature.qualifiers 	= { }
			feature.type 		= "mRNA"
			feature.id 			= geneName 
			if not geneName in genes:
				gene = SeqRecord( record.seq, annotations = {"Locus" : record.id} )
				gene.features = [ feature ]
				genes[geneName] = gene
				genes[geneName].annotations["minTarget"] = int(target.split()[1])
				genes[geneName].annotations["maxTarget"] = int(target.split()[2])
			else:
				genes[geneName].features[0].location += feature.location
				genes[geneName].annotations["minTarget"] = min([ genes[geneName].annotations["minTarget"], int(target.split()[1]) ])
				genes[geneName].annotations["maxTarget"] = max([ genes[geneName].annotations["maxTarget"], int(target.split()[2]) ])

	transMapFile.close()

	for line in transFile:
		if line.startswith('>'):
			geneName = line.split()[0][1:]
			if geneName in genes:
				length   = int(line.split()[1].split('=')[1])
				genes[geneName].annotations["startOffset"] = genes[geneName].annotations["minTarget"] - 1
				genes[geneName].annotations["endOffset"]   = length - genes[geneName].annotations["maxTarget"]

	transFile.close()

	for line in proteinFile:
		if line.startswith('>'):
			substring = line.split()[-1]
			geneName = substring.split(':')[0]
			pepStart = int(substring.split(':')[1].split('-')[0])
			pepEnd   = int(substring.split(':')[1].split('-')[1].split('(')[0] )
			pepStrand = substring.split('(')[1][0]
			print 'pepStrand', pepStrand
			
			if pepStrand == '-':
				pepStrand = -1
			else:
				pepStrand = 1

			if geneName in genes:
				geneStrand = gene.features[0].location.strand
				gene = genes[geneName]
				exons = gene.features[0]

				rcm = RangeCoordinateMapper(exons, len(gene), gene.annotations["startOffset"], gene.annotations["endOffset"] )
				
				location = rcm.rc2g(pepStart, pepEnd, pepStrand)	
				cdsFeature = SeqFeature(type = 'CDS', location = location )

				gene.features.append(cdsFeature)

	proteinFile.close()

	return genes

def 

if __name__ == '__main__':
	genes = prepareLibrary('/Users/md/marpodb/partsdb/data/map.gff3', '/Users/md/marpodb/partsdb/data/trans.fa', '/Users/md/marpodb/partsdb/data/pep.fa',  '/Users/md/marpodb/partsdb/data/genome.fa')

	for geneName, gene in genes.iteritems():
		print geneName, gene.annotations["Locus"]
		print 'SEQ', gene.seq
		for feature in gene.features:
			print feature
			print feature.extract(gene.seq)




