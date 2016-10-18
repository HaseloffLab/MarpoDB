import sys

badNames = ('Putative', 'putative', 'Predicted protein', 'Deleted.', 'Merged into', 'Uncharacterized', 'Uncharacterised')

blastpFile = open(sys.argv[1])
mapFile = open(sys.argv[2])

accList = {}

for line in mapFile:
	tabs = line.rstrip().split('\t')
	origin = tabs[5].split('(')[0].rstrip()

	if tabs[2]:
		geneName = tabs[2]
	elif tabs[3]:
		geneName = tabs[3].split()[0]
	elif tabs[4]:
		geneName = tabs[4].split()[0]
	else:
		geneName = 'N/A'

	proteinName = tabs[1].split('(')[0].rstrip()

	good = True

	for badName in badNames:
		if badName in proteinName:
			good = False
			break
	if good:
		accList[tabs[0]] = [ proteinName, origin, geneName ]

for line in blastpFile:
	tabs = line.rstrip().split('\t')
	acc = tabs[1].split('|')[1]
	if acc in accList:
		print '\t'.join( tabs + accList[acc] )

