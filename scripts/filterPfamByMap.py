import sys

inFile = open(sys.argv[1])
mapFile = open( sys.argv[2] )

mapped = []

for line in mapFile:
	if line.startswith('### '):
		transName = line.split()[1]
		if not transName in mapped:
			mapped.append(transName)


for line in inFile:
	if not line.startswith('#'):
		tabs = line.rstrip('\n').split()
		transName = tabs[3].split('|')[0]
		if transName in mapped:
			print "\t".join(tabs)
