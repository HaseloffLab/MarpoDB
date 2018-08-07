import sys

inFile = open(sys.argv[1])

for line in inFile:
	if line.startswith('##FASTA'):
		break
	if not line.startswith('#'):
		tabs = line.rstrip().split('\t')
		if tabs[2] == 'mRNA':
			keyvals = tabs[8].split(';')
			properties = { keyval.split('=')[0] : keyval.split('=')[1] for keyval in keyvals  }
			print ">{0}".format(properties["Parent"])
			print properties["translation"]

