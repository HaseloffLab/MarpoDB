import sys
from partsdb.partsdb import PartsDB
from partsdb.tools.Populators import PlantPopulator
from tables import *
from partsdb.tools.Exporters import GenBankExporter

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature

genomeFileName = 'data/genome.fa'
transcriptFileName = 'data/mapped/trans.fa'
proteinFileName = 'data/mapped/pep.fa'
mapFileName = 'data/mapped/map.gff3'


marpodb = PartsDB('postgresql:///'+sys.argv[1], clean = True, Base = Base)
marpodb.setup(prefix = "mpdb")

ppl = PlantPopulator(marpodb)

ppl.populate(mapFileName, transcriptFileName, proteinFileName,  genomeFileName)
