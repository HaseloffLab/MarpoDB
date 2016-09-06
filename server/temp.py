from marpodb import parseBlastResult
import json

inFile = open("../temp.txt")
data = inFile.read()

# j = json.loads(data)

# for key in j["BlastOutput2"][0]["report"]:
# 	print key

for row in parseBlastResult(data):
	print row