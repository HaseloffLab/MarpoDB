from Bio import SeqIO
from Bio.Seq import Seq
from compRev import compRev

def prepareLibrary(genomeFileName, transcriptFileName, proteinFileName, mapFileName):

	# Reading map file

	mapFile = open(mapFileName)

	transcripts = {}
	transcriptName = ''
	exons = []
	exonsTranscript = []
	startTarget = []
	stopTarget = []
	for line in mapFile:
		if not line.startswith('#'):
			tabs = line.split('\t')
			fType = tabs[2]
			if fType == 'cDNA_match':
				transcriptName = tabs[8].split()[0].split('=')[1]
				start = int(tabs[3])
				stop = int(tabs[4])
				direction = tabs[6]
				scaff = tabs[0]
				startTarget.append(int( tabs[8].split()[1] ) )
				stopTarget.append( int( tabs[8].split()[2] ) )
				exons.append((start,stop))
				exonsTranscript.append( (int( tabs[8].split()[1] ), int( tabs[8].split()[2] ) )  )
		
		elif line =='###\n' and transcriptName:
			start = min( [e[0] for e in exons] )
			stop  = max( [e[1] for e in exons] )
			locusName = '{0}:{1}:{2}'.format( scaff, start, stop )
		

			transcripts[transcriptName] = {}
			transcripts[transcriptName]["length"] = 0
			transcripts[transcriptName]["exons"] = []
			transcripts[transcriptName]["cdss"] = {}
			transcripts[transcriptName]["exons"] = exons
			transcripts[transcriptName]["exonsTranscript"] = exonsTranscript
			transcripts[transcriptName]["stop"] = max(stopTarget)
			transcripts[transcriptName]["start"] = min(startTarget)
			transcripts[transcriptName]["direction"] = direction
			transcripts[transcriptName]["locusName"] = locusName
			exons = []
			exonsTranscript = []
			startTarget = []
			stopTarget = []

	start = min( [e[0] for e in exons] )
	stop  = max( [e[1] for e in exons] )
	locusName = '{0}:{1}:{2}'.format( scaff, start, stop )


	transcripts[transcriptName] = {}
	transcripts[transcriptName]["length"] = 0
	transcripts[transcriptName]["exons"] = []
	transcripts[transcriptName]["cdss"] = {}
	transcripts[transcriptName]["exons"] = exons
	transcripts[transcriptName]["exonsTranscript"] = exonsTranscript
	transcripts[transcriptName]["stop"] = max(stopTarget)
	transcripts[transcriptName]["start"] = min(startTarget)
	transcripts[transcriptName]["direction"] = direction
	transcripts[transcriptName]["locusName"] = locusName

	exons = []
	startTarget = []
	stopTarget = []

	mapFile.close()

	# Get transcript lengths and sequence

	transFile = open(transcriptFileName)

	transSeq = SeqIO.to_dict(SeqIO.parse(transFile, "fasta"))

	for transcriptName, seq in transSeq.iteritems():
		if transcriptName in transcripts:
			transcripts[transcriptName]["length"] = len(seq.seq)
			# print str(seq.seq)
			# transcripts[transcriptName]["seq"] = compRev(str(seq.seq), transcripts[transcriptName]["direction"] )
			transcripts[transcriptName]["seq"] = compRev(str(seq.seq), '+' )
	transFile.close()

	# Loading scaffold sequences

	scaffFile = open(genomeFileName)

	scaffSeq = SeqIO.to_dict(SeqIO.parse(scaffFile, "fasta"))

	scaffFile.close()

	# Loading cds locations

	pepFile = open(proteinFileName)

	for line in pepFile:
		if line.startswith('>'):

			tabs = line.rstrip().split()
			transcriptName = tabs[0].split('|')[0][1:]
			cdsName 	   = tabs[0].split('|')[1]
			info           = tabs[-1]

			start = int(info.split(':')[1].split('-')[0])
			stop = int(info.split(':')[1].split('-')[1].split('(')[0])
			direction = info.split('(')[1][0]
			cdsType = tabs[5].split(':')[1]

			if transcriptName in transcripts:
				transcripts[transcriptName]["cdss"][cdsName] = {}
				transcripts[transcriptName]["cdss"][cdsName]["transLoc"] = (start, stop, direction)
				transcripts[transcriptName]["cdss"][cdsName]["geneLoc"] = []
				transcripts[transcriptName]["cdss"][cdsName]["type"] = cdsType

	pepFile.close()

	# Managing sequences

	loci = {}

	for transcriptName, transcript in transcripts.iteritems():
		locusName 	= transcript["locusName"]
		
		# Adding loci

		if not locusName in loci:

			loci[locusName] = {}

			scaff    	= locusName.split(':')[0]

			start = min( [ex[0] for ex in transcript["exons"]] )
			stop  = max( [ex[1] for ex in transcript["exons"]] )

			startCut 	= max(0, start-3000 ) + 1
			stopCut 	= min(stop + 3000, len( scaffSeq[scaff].seq ) ) + 1

			loci[locusName]["seq"] = str(scaffSeq[scaff].seq[startCut-1:stopCut]).upper()
			loci[locusName]["coordinates"] = "{0}:{1}:{2}".format(locusName.split(':')[0], startCut, stopCut)
			loci[locusName]["startCut"] = startCut
			loci[locusName]["stopCut"] = stopCut


		seq = loci[locusName]["seq"]
		startCut = loci[locusName]["startCut"]
		stopCut = loci[locusName]["stopCut"]


		# Rewriting exons position relative to a locus sequence

		for i in range( len(transcript["exons"]) ):
				exon = transcript["exons"][i]
				transcript["exons"][i] = (exon[0] - startCut + 1, exon[1] - startCut + 1 )

		
		shift = 0
		newExons = []

		print 'Old exons', transcriptName, transcript["exons"], transcript["exonsTranscript"]
		print 'Direction', transcript["direction"]
		if transcript["direction"] == '+':
			for exonG, exonT in zip(transcript["exons"], transcript["exonsTranscript"] ):
				deltaShift = (exonT[1] - exonT[0]) - (exonG[1] - exonG[0])

				print '\tGenome exon:', exonG
				print '\tTrans exon:', exonT
				print '\tShift', shift	
				seq = "".join( [ seq[ : exonG[0] - 1 + shift ], transcript["seq"][ exonT[0] - 1 : exonT[1] ], seq[ exonG[1] + shift : ]  ] )


				eGS = exonG[0] + shift
				eGE = exonG[0] + shift + (exonT[1] - exonT[0])

				newExons.append( (eGS, eGE) )

				shift += deltaShift
			transcript["exons"] = newExons
		else:
			for exonG, exonT in zip(transcript["exons"], transcript["exonsTranscript"] ):

				print '\tGenome exon:', exonG
				print '\tTrans exon:', exonT

			  	l = len(seq)

				seq = "".join( [ seq[ : exonG[0] - 1 ], compRev(transcript["seq"][ exonT[0] - 1 : exonT[1]  ], '-') , seq[ exonG[1] : ] ] )

				# print exonG[0], len(seq)

				shift = (exonT[1] - exonT[0]) - (exonG[1] - exonG[0])
				print '\tShift', shift 
				eGS = l - exonG[0] + shift
				eGE = l - exonG[1]

				newExons.append( (eGS, eGE) )

			transcript["exons"] = [ (len(seq) - e[0], len(seq) - e[1] ) for e in newExons ]
		
		transcript["sequence"] = seq
		# modifier = 0

		# for exonG, exonT in zip(transcript["exons"], transcript["exonsTranscript"] ):

		# 	newModifier = (exonT[1] - exonT[0]) - (exonG[1] - exonG[0])

		# 	seq = "{0}{1}{2}".format( seq[ :exonG[0] ] , transcript["seq"][exonT[0], exonT[1]], seq[ exonG[1]: ] )
			
		# 	exonG[0] += modifier
		# 	exonG[1]  = exonG[1] + modifier + ( exonT[1] - exonT[0] )

		# 	modifier += newModifier

	# Putting all toogether

	for transcriptName, transcript in transcripts.iteritems():
		# locusName 	= transcript["locusName"]
		
		# if not locusName in loci:

		# 	loci[locusName] = {}

		# 	scaff    	= locusName.split(':')[0]

		# 	start = min( [ex[0] for ex in transcript["exons"]] )
		# 	stop  = max( [ex[1] for ex in transcript["exons"]] )

		# 	startCut 	= max(0, start-3000 ) + 1
		# 	stopCut 	= min(stop + 3000, len( scaffSeq[scaff].seq ) ) + 1

		# 	loci[locusName]["seq"] = str(scaffSeq[scaff].seq[startCut-1:stopCut]).upper()
		# 	loci[locusName]["coordinates"] = "{0}:{1}:{2}".format(locusName.split(':')[0], startCut, stopCut)
		# 	loci[locusName]["startCut"] = startCut
		# 	loci[locusName]["stopCut"] = stopCut


		# startCut = loci[locusName]["startCut"]
		# stopCut = loci[locusName]["stopCut"]

		# for i in range( len(transcript["exons"]) ):
		# 		exon = transcript["exons"][i]
		# 		transcript["exons"][i] = (exon[0] - startCut + 1, exon[1] - startCut + 1 )
		
		for cdsName, cds in transcript["cdss"].iteritems():
				cdsDir = cds["transLoc"][2]

				started = False 
					
				if cdsDir == '+':
					cdsStart 	=  cds["transLoc"][0]
					cdsStop  	=  cds["transLoc"][1]
					cPos 		=  transcript["start"] - 1
				else:
					cdsStart 	= transcript["length"] - cds["transLoc"][1] + 1
					cdsStop  	= transcript["length"] - cds["transLoc"][0] + 1
					cPos 	 	= transcript["length"] - transcript["stop"]

				if cPos < cdsStart:

					if cdsDir == transcript["direction"]:
						transcript["exons"].sort(key = lambda x: x[0])

						for exon in transcript["exons"]:
							exonStart = exon[0]
							exonStop  = exon[1]
					
							exonL = exonStop - exonStart + 1
							cPos += exonL

							if (not started):
								if (cPos >= cdsStart):
									if (cPos < cdsStop):
										cds["geneLoc"].append( ( exonStop - (cPos - cdsStart), exonStop ) )
										started = True
									else:
										cds["geneLoc"].append( ( exonStop - (cPos - cdsStart), exonStop - (cPos - cdsStop) ) )
										break
							else:
								if (cPos < cdsStop):
									cds["geneLoc"].append( ( exonStart, exonStop ) )
								else:
									cds["geneLoc"].append( ( exonStart, exonStop - (cPos - cdsStop) ) )
									break
					else:
						transcript["exons"].sort(key = lambda x: x[0], reverse = True)

						for exon in transcript["exons"]:
							exonStart = exon[0]
							exonStop  = exon[1]
					
							exonL = exonStop - exonStart + 1
							cPos += exonL

							if (not started):
								if (cPos >= cdsStart):
									if (cPos < cdsStop):
										cds["geneLoc"].append( ( exonStart, exonStart + (cPos - cdsStart) ) )
										started = True
									else:
										cds["geneLoc"].append( ( exonStart + (cPos - cdsStop), exonStart + (cPos - cdsStart) ) )
										break
							else:
								if (cPos < cdsStop):
									cds["geneLoc"].append( ( exonStart, exonStop ) )
								else:
									cds["geneLoc"].append( ( exonStart + (cPos - cdsStop), exonStop ) )
									break

	return loci, transcripts

	# conn = psycopg2.connect("dbname={0}".format(dbName))
	# cur = conn.cursor()

	# commonPrefix = 'mpdb'

	# def newID(table):
	# 	cur.execute("SELECT COUNT(*) FROM {0}".format(table), )
	# 	nid = cur.fetchone()[0]
	# 	if nid is not None:
	# 		nid = "{0}.{1}.{2}".format(commonPrefix, table, int(nid) )
	# 		return nid
	# 	else:
	# 		return None


	# def insertRow(table, valuesDict):
	# 	nid = newID(table)

	# 	if nid:
	# 		valuesDict['id'] = nid
	# 		keys = tuple( valuesDict.keys() )
	# 		values = tuple( [ valuesDict[key] for key in keys ] )
			
	# 		keyString = ','.join(keys)

	# 		exeString = "INSERT INTO {0} ({1}) VALUES %s RETURNING id".format(table, keyString)
	# 		cur.execute(exeString, (values,))
			

	# 		cid = cur.fetchone()[0]

	# 		if cid:
	# 			conn.commit()
	# 			return cid

	# 	return None

	# def compRev(seq, sign):
	# 	if sign =='-':
	# 		sub = dict( zip( ['A', 'T', 'G', 'C', 'N'], ['T', 'A', 'C', 'G', 'N'] ) )
	# 		return ''.join( sub[nucl] for nucl in seq.upper() )[::-1]
	# 	else:
	# 		return seq

	# def getKey(dic, value):
	# 	ret = None
	# 	for cid, cval in dic.iteritems():
	# 		if cval == value:
	# 			ret = cid
	# 	return ret

	# def invertDir(d):
	# 	if d == '+':
	# 		return '-'
	# 	else:
	# 		return '+'

	# def getUtrCoordinates(exons, cdsStart, cdsStop):
	# 	coordinates = []

	# 	for exon in sorted(exons, key = lambda x: x[0]):
	# 		if cdsStart < exon[0]:
	# 			if cdsStop <= exon[1] and cdsStop >= exon[0]:
	# 				coordinates.append( (exon[0], cdsStop) )
	# 				break
	# 			elif cdsStop > exon[1]:
	# 				coordinates.append( exon )
	# 		elif cdsStart <= exon[1]:
	# 			if cdsStop <= exon[1]:
	# 				coordinates.append( (cdsStart, cdsStop) )
	# 				break
	# 			else:
	# 				coordinates.append( (cdsStart, exon[1]) )

	# 	return coordinates

	# for locusName, locus in loci.iteritems():
	# 	locus["PromoterP"] = None
	# 	locus["PromoterN"] = None
	# 	locus["TerminatorP"] = None
	# 	locus["TerminatorN"] = None
	# 	locus["ID"] = insertRow('locus', {'coordinates' : locus["coordinates"]})

	# for transcriptName, transcript in transcripts.iteritems():

	# 	locus = loci[ transcript["locusName"] ]

	# 	transStart = min ( [ex[0] for ex in transcript["exons"]] )
	# 	transStop =  max ( [ex[1] for ex in transcript["exons"]] )

	# 	for cdsName, cds in transcript["cdss"].iteritems():
	# 		if transcript["direction"] == cds["transLoc"][2]:
	# 			direction = '+'
	# 		else:
	# 			direction = '-'

	# 		coordinates = ';'.join(str(cds[0]) + ',' + str(cds[1]) + ',' + direction for cds in cds["geneLoc"])

	# 		cdsStart = min( [ex[0] for ex in  cds["geneLoc"]] )
	# 		cdsStop  = max( [ex[1] for ex in  cds["geneLoc"]] )


	# 		sequence = compRev(locus["seq"][cdsStart-1:cdsStop], direction)

	# 		cdsID = insertRow('cds', {'coordinates' : coordinates, 'seq' : sequence})


	# 		if direction == '+':
	# 			if not locus["PromoterP"]:
	# 				sequence = compRev(locus["seq"][:transStart-1], direction)
	# 				promoterID = insertRow('promoter', {'seq' : sequence})

	# 				sequence = compRev(locus["seq"][transStop-1:], direction)
	# 				terminatorID = insertRow('terminator', {'seq' : sequence})
	# 				locus["PromoterP"] = promoterID
	# 				locus["TerminatorP"] = terminatorID
	# 			else:
	# 				promoterID = locus["PromoterP"]
	# 				terminatorID = locus["TerminatorP"]
			
	# 		if direction == '-':
	# 			if not locus["PromoterN"]:
	# 				sequence = compRev(locus["seq"][transStop-1:], direction)
	# 				promoterID = insertRow('promoter', {'seq' : sequence})

	# 				sequence = compRev(locus["seq"][:transStart-1], direction)
	# 				terminatorID = insertRow('terminator', {'seq' : sequence})
	# 				locus["PromoterN"] = promoterID
	# 				locus["TerminatorN"] = terminatorID
	# 			else:
	# 				promoterID = locus["PromoterN"]
	# 				terminatorID = locus["TerminatorN"]

	# 		utrLD = {'+':'utr5', '-':'utr3'}
	# 		utrRD = {'+':'utr3', '-':'utr5'}

	# 		sequence = compRev(locus["seq"][transStart-1:cdsStart], direction)
	# 		coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)

	# 		utrL = insertRow(utrLD[direction] , {'seq': sequence, 'coordinates': coordinates})

	# 		sequence = compRev(locus["seq"][cdsStop-1:transStop], direction)
	# 		coordinates = getUtrCoordinates(cds["geneLoc"], cdsStart, cdsStop)

	# 		utrR = insertRow(utrRD[direction], {'seq': sequence, 'coordinates': coordinates})

	# 		geneID = insertRow('gene', {'promoter_id' : promoterID, '{0}_id'.format( utrLD[direction] ) : utrL, '{0}_id'.format( utrRD[direction] ) : utrR, 'cds_id' : cdsID, 'terminator_id' : terminatorID, "locus_id" : locus["ID"]  })

	# 		if geneID:
	# 			if direction == '+':
	# 				print '\t'.join([cdsID, geneID, promoterID, utrL, cdsID, utrR, terminatorID, locus["ID"], "{0}|{1}".format(transcriptName, cdsName) ])
	# 			else:
	# 				print '\t'.join([cdsID, geneID, promoterID, utrR, cdsID, utrL, terminatorID, locus["ID"], "{0}|{1}".format(transcriptName, cdsName) ])
