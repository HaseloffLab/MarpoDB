from partsdb.partsdb import PartsDB
from partsdb.tools.Populators import PlantPopulator
from server.tables import *
from partsdb.tools.Exporters import GenBankExporter

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature

genomeFileName = 'data/genome.fa'
transcriptFileName = 'data/trans.fa'
proteinFileName = 'data/pep.fa'
mapFileName = 'data/map.gff3'


marpodb = PartsDB('postgresql:///testdb', clean = True, Base = Base)
marpodb.setup(prefix = "mpdb")

ppl = PlantPopulator(marpodb)

ppl.populate(mapFileName, transcriptFileName, proteinFileName,  genomeFileName)
