import sys
import os

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from utils import *
from marpodb.client import MarpoDBClient
from marpodb.client.tables import *

client = MarpoDBClient("/marpodbtak")
session = client.Session()

inputList = open("../../data/synthesis/21_05_2018_TF_SYNTHESIS.tsv")
outPath = "../../output/synthesis/"

outputList = open( os.path.join(outPath, "outputList.tsv"), "w" )

def domesticate(rec, type, name, description, reference, author, group, dbid, addBsaI = False):
	# Adding overhangs
	rec = Seq(rss[parts[type][0]], IUPAC.unambiguous_dna) + rec + Seq(rss[parts[type][1]], IUPAC.unambiguous_dna)
	
	# Removing restriction sites
	rec = recode(rec, type, name)
	
	if not rec:
		return None

	# Adding annotation

	for i in range( len(rec.features) ):
		rec.features[i].qualifiers["ApEinfo_fwdcolor"] = colours[ i%10 ] 

	feature = SeqFeature( FeatureLocation(rsLen, len(rec)-rsLen),
		type = 'misc_feature', qualifiers = {'label' : labels[type] + name, "ApEinfo_fwdcolor" : partColors[type] } )
	
	over5   = SeqFeature( FeatureLocation(0, rsLen),
		type='misc_feature', id = 'Overhang',
			qualifiers = {'label' : type + ' overhang',  "ApEinfo_fwdcolor" : "#56FA76" } )
	over3   = SeqFeature( FeatureLocation(len(rec) - rsLen, len(rec) ),
		type='misc_feature', id = 'Overhang', qualifiers = {'label' : type + ' overhang', "ApEinfo_fwdcolor" : "#56FA76" } )


	rec.features.insert(0,feature)
	rec.features.append(over5)
	rec.features.append(over3)

	# Adding BsaI
	if addBsaI:
		rec = Seq("GGTCTCA", IUPAC.unambiguous_dna) + rec + Seq("CGAGACC", IUPAC.unambiguous_dna)

	recauthor = Reference()
	
	recauthor.authors = author + "({0})".format(group)
	recauthor.title = 'Marchanita polymorpha synthesis'
	recauthor.journal = 'Haseloff Lab'

	rec.annotations["references"] = [ recauthor ]

	if reference.strip():
		ref = Reference()
		ref.authors = ''
		ref.title = ''
		ref.journal = ''
		ref.comment = reference
		rec.annotations["references"].append(ref)

	rec.description = description
	rec.name = name[0:16]
	rec.annotations["accession"] = ''
	rec.annotations["version"] = ''

	if dbid.strip():
		rec.dbxrefs = [ dbid ]

	return rec

def processPart(rec, alias, type, tff, infoName):
	recd = domesticate(rec, type = type, name = alias, description = tff + " transcription factor", reference = "", author = "M. Delmans", group = "Haseloff Lab", dbid = alias, addBsaI = True)
	if recd:
		fileName = "{0}_{1}.gb".format(type, alias)
		outFile = os.path.join(outPath, fileName)
		SeqIO.write(recd, outFile, "genbank")
		return fileName
	else:
		return "Domestication failed"

outputList.write("\t".join(["alias", "tff", "INFO_NAME", "status", "PROM5", "PROM", "UTR5", "PROM5_LEN", "PROM_LEN", "5UTR_LEN"]) + "\n")

for line in inputList:
	if not line.startswith("#"):
		tabs = line.split('\t')

		alias = tabs[0]
		tff = tabs[1]
		infoName = tabs[2].rstrip()

		gene = session.query(Gene).filter(Gene.alias == alias).first()

		status = "OK"
		prom5 =""
		prom = ""
		utr5 = ""

		utr5Len		= "0"
		promLen		= "0"
		prom5Len	= "0"

		if gene:
			promoterSeq = gene.promoter.seq.upper()
			
			# Doing promoter checks
			if len(promoterSeq) < minPromoterLen:
				status = "Promoter is smaller than {0}".format(minPromoterLen)
				outputList.write("\t".join([ alias, tff, infoName, status, prom5, prom, utr5, prom5Len, promLen, utr5Len ]) + "\n")
				continue
			if promoterSeq[-minPromoterLen:].count('N'):
				status = "Promoter has Ns in first {0} bp".format(minPromoterLen)
				outputList.write("\t".join([ alias, tff, infoName, status, prom5, prom, utr5, prom5Len, promLen, utr5Len ]) + "\n")
				continue

			utr5Seq = gene.utr5.seq.upper() if gene.utr5 else ""

			if utr5Seq.count('N') > 0:
				status = "UTR5 has Ns"
				outputList.write("\t".join([ alias, tff, infoName, status, prom5, prom, utr5, prom5Len, promLen, utr5Len ]) + "\n")
				continue

			promExtract = extractPROM5(gene, alias)
			
			if len(promExtract) == 1:
				prom5 = processPart( promExtract[0], alias, "PROM5", tff, infoName )
				prom5Len = str( len(promExtract[0].seq) + 8 )

			else:
				prom = processPart( promExtract[0], alias, "PROM", tff, infoName )
				utr5 = processPart( promExtract[1], alias, "5UTR", tff, infoName )
				promLen  = str( len(promExtract[0].seq) + 8 )
				utr5Len = str( len(promExtract[1].seq) + 8 )

			outputList.write("\t".join([ alias, tff, infoName, status, prom5, prom, utr5, prom5Len, promLen, utr5Len ]) + "\n")