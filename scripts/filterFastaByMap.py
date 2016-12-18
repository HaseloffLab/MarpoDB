import sys

inFile = open(sys.argv[1])
mapFile = open( sys.argv[2] )

mapped = []

for line in mapFile:
	if line.startswith('### '):
		transName = line.split()[1].split('|')[0]
		if not transName in mapped:
			mapped.append(transName)

keep = False

for line in inFile:
	line = line.rstrip('\n')
	if line.startswith('>'):
		tabs = line.split()
		transName = tabs[0].split('|')[0][1:]
		if transName in mapped:
			keep = True
		else:
			keep = False

	if keep:
		print line
