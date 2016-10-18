# Filtering blast file (arg1) by coverage (arg2) and identity (arg3)

import sys

inFile = open(sys.argv[1])
thCov = int(sys.argv[2])
thId = float(sys.argv[3])

mapped = []

for line in inFile:
	tabs = line.rstrip('\n').split('\t')
	try:
		qcov = int(tabs[8])
		hid = float(tabs[9])

		if qcov > thCov and hid > thId:
				print '\t'.join(tabs)
	except:
		pass