from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference, CompoundLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

maxSynLen 		= 1800
minPromoterLen 	= 1000
minSynLen		= 0

rsLen = 4

maxPartLen = maxSynLen - 2 * (rsLen)

# A			  X           B      C      D      E      F      G      H      I
#	PromoterD   PromoterP   5UTR   NTAG   CDS1   CDS2   CTAG   3UTR   TERM

rss ={	'A': 'GGAG',
		'B': 'TACT',
		'C': 'CCAT',
		'D': 'AATG',
		'E': 'AGCC',
		'F': 'TTCG',
		'G': 'GCTT',
		'H': 'GGTA',
		'I': 'CGCT',
		'X': 'TGAC'}

parts = {	'PROM5' : ('A', 'D'),
			'PROM'  : ('A', 'B'),
			'5UTR'	: ('B', 'D'),
			'CDS'	: ('D', 'G'),
			'CDS12' : ('D', 'F'),
			'CDS1'	: ('D', 'E'),
			'CDS2'	: ('E', 'F'),
			'CTAG'	: ('F', 'G'),
			'TERM'  : ('H', 'I'),
			'3TERM' : ('G', 'I'),
			'TU'	: ('A', 'I'),
			'PROMD'	: ('A', 'X'),
			'PROMP' : ('X', 'B')
}

labels = {
			'PROM5' : "pro",
			'PROM'  : "pro",
			'5UTR'	: "utr5",
			'CDS'	: "",
			'CDS12'  : "",
			'CDS1'	: "",
			'CDS2'	: "",
			'-CDS'  : "",
			'CTAG'	: "ctag",
			'TERM'  : "term",
			'3TERM' : "3term",
			'TU'	: "",
			'PROMD' : "prod",
			'PROMP' : "prop"
}

def coordinatesToLocation(coordinates):
	locationParts = [ FeatureLocation(int(p[0]), int(p[1]), int(p[2]) ) for p in [ s.split(',') for s in coordinates.split(';')] ]
	if len(locationParts) == 1:
		return locationParts[0]
	elif len(locationParts) > 1:
		return CompoundLocation(locationParts)
	else:
		return None


def inFrame(seq):
	return ( len(seq) % 3 == 0, "CDS not in frame." )

def noStop(seq):
	return ( seq.translate().find('*') == -1, "CDS terminates too early.")

def ends(seq):
	return ( seq[-3:].translate() == '*', "CDS lacks stop codon." )

cdsCheckFlanks = {
					"CDS"   : ["ATG", ""],
					"CDS12" : ["ATG", "T"],
					"CDS1"	: ["ATG", "A"],
					"CDS2"  : ["GCC", "T"],
					"CTAG"	: ["TCG", ""]
				}
cdsChecks = {
				"CDS" 	: [ inFrame, lambda rec: noStop(rec[:-3]), ends ],
				"CDS12" : [ inFrame, noStop ],
				"CDS1"  : [ inFrame, noStop ],
				"CDS2"	: [ inFrame, noStop ],
				"CTAG"	: [ inFrame, lambda rec: noStop(rec[:-3]), ends ]
}

def prepareCDS(rec, type, name):
	# Validate sequence
	for check in cdsChecks[type]:
		res = check( cdsCheckFlanks[type][0] + rec.seq + cdsCheckFlanks[type][1] )
		if not res[0]:
			print "{0}_{1}:\t{2}".format(type, name, res[1])

	if len(rec) > maxPartLen:
		print '{0}_{1}:\tToo long'.format(type, name)
		return None

	return rec


def prepare5(rec, type, name):
	dif = len(rec) - maxPartLen
	if dif > 0:
		rec = rec[dif:]
		print "{0}_{1}:\tShortened by {2}nt from 5' end.".format(type,name,dif)
	return rec	

def prepare3(rec, type, name):
	dif = len(rec) - maxPartLen
	if dif > 0:
		rec = rec[:len(rec)-dif]
		print "{0}_{1}:\tShortened by {2}nt from 3' end.".format(type,name,dif)
	return rec

def prepareM(rec, type, name):
	if len(rec) > maxPartLen:
		print '{0}_{1}:\tToo long'.format(type, name)
		return None
	return rec

prepareFuns = {
			'PROM5' : prepare5,
			'PROM'  : prepare5,
			'PROMD' : prepare5,
			'PROMP' : prepare5,
			'5UTR'	: prepareM,
			'CDS'	: prepareCDS,
			'CDS12' : prepareCDS,
			'CDS1'	: prepareCDS,
			'CDS2'	: prepareCDS,
			'CTAG'	: prepareCDS,
			'TERM'  : prepare3,
			'3TERM' : prepare3,
			'TU'	: prepare5
}

extract = {
			'PROM5' : lambda gene: Seq(gene.promoter.seq + (gene.utr5.seq if gene.utr5 else ''), IUPAC.unambiguous_dna),
			'PROM'  : lambda gene: Seq(gene.promoter.seq, IUPAC.unambiguous_dna),
			'5UTR'	: lambda gene: Seq(gene.utr5.seq, IUPAC.unambiguous_dna) if gene.utr5 else '',
			'CDS'	: lambda gene: coordinatesToLocation(gene.cds.coordinates).extract( Seq(gene.cds.seq, IUPAC.unambiguous_dna)[3:] ),
			'CDS12' : lambda gene: coordinatesToLocation(gene.cds.coordinates).extract( Seq(gene.cds.seq, IUPAC.unambiguous_dna)[3:-3] + "GC" ),
			'CDS1'	: lambda gene: coordinatesToLocation(gene.cds.coordinates).extract( Seq(gene.cds.seq, IUPAC.unambiguous_dna)[3:-3] + "GC" ),
			'CDS2'	: lambda gene: coordinatesToLocation(gene.cds.coordinates).extract( Seq(gene.cds.seq, IUPAC.unambiguous_dna)[3:-3] + "GC" ),
			'CTAG'	: lambda gene: coordinatesToLocation(gene.cds.coordinates).extract( Seq(gene.cds.seq, IUPAC.unambiguous_dna)[3:] ),
			'TERM'  : lambda gene: Seq(gene.terminator.seq, IUPAC.unambiguous_dna),
			'3TERM' : lambda gene: Seq((gene.utr3.seq if gene.utr3 else '') + gene.terminator.seq, IUPAC.unambiguous_dna),
			'TU'	: lambda gene: None,
			'PROMD' : lambda gene: None,
			'PROMP' : lambda gene: None
}

## Adaptation of synthesis_recode.py (Bernardo Pollak)
colours =['#99ccff','#F96381','#ffff80','#ff9966','#ccb3ff','#4db8ff','#ffcc66','#d5ff80','#C6A49A','#56FAAF']

replace = {
			"GGTCTC" 	: ["GGaCTC","GGTCcC","GGTaTC"],
			"GAGACC" 	: ["GAaACC","GAGAtC" ,"GAGgCC"],
			"GCTCTTC"	: ["GCcCTTC","GCTCgTC","GCTgTTC"],
			"GAAGAGC"	: ["GAgGAGC","GAAGgGC","GAAaAGC"]
}

def recfind(pattern, string, where_should_I_start=0):
    # Save the result in a variable to avoid doing the same thing twice
    pos = string.find(pattern, where_should_I_start)
    if pos == -1:
        # Not found!
        return []
    # No need for an else statement
    return [pos] + recfind(pattern, string, pos + len(pattern))

def check(seq, sites):
	found = 0
	for site in sites:
		a = seq.find(site, 0, len(seq))
		if a != -1:
			found = found + 1
	return found

def recode(rec, type, name):

	sites = {}

	# Finding restriction sites
	for pattern in replace:
		sites[pattern] = recfind(pattern, rec.seq.upper())

	# Replacing and labeling restriction sites
	n = 0
	for pattern in sites:
		for site in sites[pattern]:
			rec = rec[:site] + Seq( replace[pattern][ (site-1) % 3 ], IUPAC.unambiguous_dna ) + rec[site + len(pattern):]
			
			i = n % 10
			n = n + 1

			recFeature = SeqFeature( FeatureLocation(site, site + len(pattern)),
				type = "misc_feature", id = "Restriction site",
					qualifiers = {"label" : "Restriction site removed", "ApEinfo_fwdcolor": colours[i] } )
			rec.features.append(recFeature)

	# Checking there are no more sites present
	if check(rec.seq.upper(), sites) > 0:
		print "{0}_{1}:\tRestriction sites still present after replacement.".format(type, name)
		return None

	return rec

partColors = {
	'PROM5' : "#ffff4d",
	'PROM'  : "#ffff4d",
	'5UTR'	: "#5073ff",
	'CDS'	: "#3385ff",
	'CDS12' : "#3385ff",
	'CDS1'	: "#3385ff",
	'CDS2'	: "#3385ff",
	'CTAG'	: "#3385ff",
	'TERM'  : "#ff5050",
	'3TERM' : "#ff5050",
	'TU'	: "#5cd65c",
	'PROMD' : "#ffff4d",
	'PROMP' : "#ffff4d"
}

## END of Adaptation of synthesis_recode (Bernardo Pollak)

def domesticate(rec, type, name, description, reference, author, group, dbid, addBsaI = False):

	# Cutting sequences to fit the size, checking CDSs
	rec =  prepareFuns[type](rec, type, name)
	
	if not rec:
		return None
	
	# Adding overhangs
	rec = Seq(rss[parts[type][0]], IUPAC.unambiguous_dna) + rec + Seq(rss[parts[type][1]], IUPAC.unambiguous_dna)
	
	# Removing restriction sites
	rec = recode(rec, type, name)
	
	if not rec:
		return None

	if len(rec) < minSynLen:
		print "{0}_{1}:\tSequence too short.".format(type, name)
		return None

	rec.seq = rec.seq.upper()

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

def extractPROM5(gene, name):
	if gene.utr5:
		if len(gene.utr5.seq) <= maxPartLen - minPromoterLen:
			seq = gene.promoter.seq + gene.utr5.seq
			cutStart = max( len(seq) - maxPartLen, 0 )
			seq = seq[cutStart:]

			rec = SeqRecord(seq = Seq(seq.upper(), IUPAC.unambiguous_dna), name = name)

			promEnd = len(gene.promoter.seq) - cutStart

			rec.features.append( SeqFeature( FeatureLocation(0, promEnd)		, type = 'PROM', qualifiers = {'label' : 'pro'  + name }))
			rec.features.append( SeqFeature( FeatureLocation(promEnd, len(seq))	, type = '5UTR', qualifiers = {'label' : 'utr5' + name }))
			
			return [ rec ]

		else:
			cutStart = max( len(gene.promoter.seq) - maxPartLen, 0 )
			promoterSeq = gene.promoter.seq[ cutStart: ]
			promoterRec = SeqRecord(seq = Seq(promoterSeq.upper(), IUPAC.unambiguous_dna), name = name)

			utrRec = SeqRecord(seq = Seq( gene.utr5.seq.upper(), IUPAC.unambiguous_dna), name = name)

			return [ promoterRec, utrRec ]

	else:
		cutStart = max( len(gene.promoter.seq) - maxPartLen, 0 )
		seq = gene.promoter.seq[cutStart:]
		rec = SeqRecord(seq = Seq(seq.upper(), IUPAC.unambiguous_dna), name = name)

		return [ rec ]