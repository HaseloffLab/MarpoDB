import sys
from partsdb.partsdb import PartsDB
from tables import *

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

marpodb.annotate('blastphit', fileName = 'data/filtered/blastp.info')
marpodb.annotate('pfamhit', fileName = 'data/filtered/Pfam.domtblout')
