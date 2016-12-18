 Filtering blast file (arg1) by identity (arg2)

import sys

inFile = open(sys.argv[1])
thCov = float(sys.argv[2])
thId = float(sys.argv[3])


for line in inFile:
	tabs = line.rstrip('\n').split('\t')

	cov = 100*(float(tabs[5])-float(tabs[4])) / float(tabs[2])
	hid = float(tabs[8])
	if hid > thId and cov > thCov:
		print '\t'.join(tabs[0:8])+'\t'+str(cov)+'\t'+'\t'.join(tabs[8:])
#		print '\t'.join(tabs)