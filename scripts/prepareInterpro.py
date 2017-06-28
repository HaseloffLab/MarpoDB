import sys
from BCBio import GFF

listFile = open(sys.argv[1])

description = {}

for line in listFile:
	tabs = line.rstrip().split()
	if tabs[0].startswith("IPR"):
		description[tabs[0]] = " ".join(tabs[1:])

gffFile = open(sys.argv[2])

headers = ["dbid", "origin", "start", "end", "description", "refID", "eVal"]

print "\t".join(headers)

for line in gffFile:
	if line.startswith("##FASTA"):
		break
	if not line.startswith('#'):
		tabs 	= line.split('\t')
		mpdb 	= tabs[0]
		type 	= tabs[1]
		start 	= int(tabs[3])
		end 	= int(tabs[4])

		props = tabs[8].split(";")

		try:
			keys = [ prop.split("=")[0] for prop in props ]
			vals = [ prop.split("=")[1].split(",") for prop in props ]
		except:
			print "Error parsing"

		desc = { key:val for key, val in zip( keys, vals ) }
		if "Dbxref" in desc:
			
			eVal = 1			
			refId = desc["Dbxref"][0].split(":")[1].strip()
			refId = refId.replace('"', '')

			if not refId in description:
				continue

			if tabs[5] != ".":
				eVal = float(tabs[5])
				if eVal > 1:
					eVal = 1

			if "Ontology_term" in desc:
				for goTerm in desc["Ontology_term"]:
					print "\t".join([mpdb, "GO", "0", "0", goTerm[1:-1], goTerm[1:-1], "1" ])

			print "\t".join([mpdb, type, str(start), str(end), description[refId], refId, str(eVal)])
			# break