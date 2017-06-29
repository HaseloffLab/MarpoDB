import sys
from partsdb.partsdb import PartsDB
from tables import *

marpodb = PartsDB('postgresql:///'+sys.argv[1], Base = Base)

# marpodb.annotate('blastphit', fileName = 'data/filtered/blastp.info.small')
# marpodb.annotate('interprohit', fileName = 'data/interpro.list.small')
marpodb.annotate('dbxref', fileName = 'data/DbxRef.list')