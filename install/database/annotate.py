import sys
from partsdb.partsdb import PartsDB
from server.tables import *

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

marpodb.annotate('blastphit', fileName = 'data/filtered/blastp_filtered.outfmt6')
marpodb.annotate('pfamhit', fileName = 'data/filtered/Pfam.domtblout')
