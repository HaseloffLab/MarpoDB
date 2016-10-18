import sys

for line in open(sys.argv[1]):
	line = line.rstrip()
	if line.startswith('>'):
		line = line.replace('|', ' ')
		line = line[0] + 'lcl|' + line[1:]
	print line
