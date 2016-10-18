import sys

transcripts = []

mRNAFile = open(sys.argv[1] + '.transdecoder.mRNA')
mRNAOutFile = open(sys.argv[1] + '.transdecoder.complete.mRNA', 'w')

for line in mRNAFile:
	if line.startswith('>'):
		tabs = line.split()
		if tabs[5] == "type:complete":
			keep = True
			
			transName = tabs[0].split('|')[0][1:]
			if not transName in transcripts:
				transcripts.append(transName)

		else:
			keep = False

	if keep:
		mRNAOutFile.write(line)

mRNAFile.close()
mRNAOutFile.close()

transcriptFile = open(sys.argv[2])
transcriptOutFile = open(sys.argv[1] + '.transdecoder.complete.trans', 'w')

for line in transcriptFile:
	if line.startswith('>'):
		transName = line.split()[0][1:]
		if transName in transcripts:
			keep = True
		else:
			keep = False
	if keep:
		transcriptOutFile.write(line)

transcriptOutFile.close()


pepFile = open(sys.argv[1] + '.transdecoder.pep')
pepOutFile = open(sys.argv[1] + '.transdecoder.complete.pep', 'w')

for line in pepFile:
	if line.startswith('>'):
		tabs = line.split()
		if tabs[5] == "type:complete":
			keep = True
		else:
			keep = False
	if keep:
		pepOutFile.write(line)

pepFile.close()
pepOutFile.close()


