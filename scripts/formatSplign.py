import sys
import re

queries = {}

def sign(x):
	if x > 0:
		return '+'
	else:
		return '-'

exonPattern = re.compile(".*<exon>.*")

for line in open(sys.argv[1]):

	if not line.startswith('#'):

		tabs = line.rstrip().split()
		qId = abs( int(tabs[0]) )
		qStrand = sign( int(tabs[0]) )

		if not qId in queries:
			queries[qId] = {}
			queries[qId]['qStrand'] = qStrand
			queries[qId]['exons'] = []
			queries[qId]['length'] = 0
			queries[qId]['covered'] = 0
			queries[qId]['identified'] = 0
			queries[qId]['scaff'] = tabs[2]
			queries[qId]['target'] = tabs[1]
			queries[qId]['code'] = ''

		if queries[qId]['qStrand'] == qStrand:
			length = int(tabs[4])
			queries[qId]['length'] += length

			if exonPattern.match(tabs[9]):
				queries[qId]['code'] += tabs[10]
				queries[qId]['covered'] += length
				identity = float(tabs[3])
				queries[qId]['identified'] += length * identity

				tStart = int(tabs[7])
				tStop = int(tabs[8])

				if tStart < tStop:
					queries[qId]['tStrand'] = '+'
				else:
					queries[qId]['tStrand'] = '-'

				exonQ = sorted ( ( int(tabs[5]), int(tabs[6]) ), key = lambda x: x )
				exonT = sorted ( ( int(tabs[7]), int(tabs[8]) ), key = lambda x: x )
				queries[qId]['exons'].append( (exonQ, exonT) )


coverageTh = float(sys.argv[2])
identityTh = float(sys.argv[3])

for qId in queries:
	coverage = queries[qId]['covered'] / float(queries[qId]['length'])
	identity = queries[qId]['identified'] / float(queries[qId]['length'])

	if (coverage > coverageTh) and (identity > identityTh) and (not 'D' in queries[qId]['code']) and (not 'I' in queries[qId]['code']):
		print '###'
		print '###', queries[qId]['target'], coverage, identity, queries[qId]['code']
		queries[qId]['exons'].sort( key = lambda x: x[0][0] )
		for exon in queries[qId]['exons']:
			tabs = []
			tabs.append(queries[qId]['scaff'])
			tabs.append('Cam1')
			tabs.append('cDNA_match')
			tabs.append(str(exon[1][0]) )
			tabs.append(str(exon[1][1]) )
			tabs.append('.')
			
			if queries[qId]['qStrand'] == queries[qId]['tStrand']:
				tabs.append('+')
			else:
				tabs.append('-')

			tabs.append('.')
			tabs.append("Target={0} {1} {2}".format(queries[qId]['target'], exon[0][0], exon[0][1]))
			print '\t'.join(tabs)	
		
	