 # Populate marpodb from tak genome files
import sys

from marpodb.core import MarpoDB
from marpodb.core.populators import TakPopulator

marpodb = MarpoDB("/marpodb4")
marpodb.setup(prefix="mpdb")

populator = TakPopulator(marpodb)

if len(sys.argv) == 3:
	genomeFileName = sys.argv[1]
	mapFileName = sys.argv[2]
else:
	genomeFileName = "../../data/tak3.1/Mpolymorpha_320_v3.0.fa"
	mapFileName = "../../data/tak3.1/Mpolymorpha_320_v3.1.gene.gff3"

print "Parsing genes..."

populator.populate(mapFileName, genomeFileName, datasetName = "tak", datasetVersion = "3.1")