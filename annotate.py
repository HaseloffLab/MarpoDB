from partsdb.partsdb import PartsDB
from server.tables import *

marpodb = PartsDB('postgresql:///testdb', Base = Base)

marpodb.annotate('blastphit', fileName = 'data/blastp_sample.info')
marpodb.annotate('pfamhit', fileName = 'data/Pfam_sample.domtblout')