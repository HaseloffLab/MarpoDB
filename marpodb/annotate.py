from partsdb.partsdb import PartsDB
from server.tables import *

marpodb = PartsDB('postgresql:///testdb', Base = Base)

marpodb.annotate('blastphit', fileName = 'data/marpodb_pep.info')
marpodb.annotate('pfamhit', fileName = 'data/marpodb_pep.pfam.domtblout')