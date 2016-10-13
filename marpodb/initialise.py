from partsdb.partsdb import PartsDB
from partsdb.tools.Populators import PlantPopulator
from tables import *
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
# exporter = GenBankExporter(marpodb)

# for gene in marpodb.session.query(Gene).all():
# 	feature = SeqFeature( location = exporter.coordinatesToLocation(gene.cds.coordinates) )
# 	seq = feature.extract(Seq(gene.cds.seq, generic_dna)).translate()
# 	record = SeqRecord( seq = seq, id = gene.cds.dbid )
# 	records.append(record)

# outputFile = open('CDSs.fa', 'w')

# SeqIO.write(records, outputFile, 'fasta')