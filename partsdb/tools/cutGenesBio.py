from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from pprint import pprint
from Bio.SeqFeature import SeqFeature


def prepareLibrary(transMapFileName, cdsMapFileName, genomeFile):
	
	transMapFile = open(transMapFileName)
	cdsMapFile   = open(cdsMapFileName)
	genomeFile = open(genomeFile)
	
	scaffoldDict = SeqIO.to_dict( SeqIO.parse(genomeFile, "fasta") )
	genomeFile.close()

	transMapRecords = GFF.parse(transMapFile, base_dict = scaffoldDict)
	cdsMapRecords = GFF.parse(cdsMapFile)

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

			genes[geneName].features[0].location += feature.location
	transMapFile.close()

	for record in cdsMapRecords:
		for feature in record.features:
			target = feature.qualifiers["Target"][0]
			geneName  = '_'.join(target.split()[0].split('_')[0:2])

			feature.qualifiers 	= { }
			feature.type 		= "CDS"
			feature.id 			= target.split()[0]

			if geneName in genes:
				gene = genes[geneName]

				found = False
				for gfeature in gene.features:
					if gfeature.id == feature.id:
						gfeature.location += feature.location
						found = True
						break
				
				if not found:
					gene.features.append(feature)
	
	cdsMapFile.close()



	loci = {}


	for geneName, gene in genes.iteritems():

		start = min( [feature.location.start 	for feature in gene.features] )
		end   = max( [feature.location.end 		for feature in gene.features] )
		
		cutStart 	= max( [start-3000, 0] )
		cutEnd		= min( [end+3000, len(gene.seq)] )

		locusName = "{0}:{1}-{2}".format(gene.annotations["Locus"], cutStart, cutEnd )
		genes[geneName] = gene[cutStart:cutEnd]
		genes[geneName].annotations["Locus"] = locusName

	return genes

if __name__ == '__main__':
	genes = prepareLibrary('/Users/md/marpodb/partsdb/data/map.gff3', '/Users/md/marpodb/partsdb/data/cds.gff3',  '/Users/md/marpodb/partsdb/data/genome.fa')

	for geneName, gene in genes.iteritems():
		print geneName, gene.annotations["Locus"]
		for feature in gene.features:
			print feature




