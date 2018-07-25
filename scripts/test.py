from marpodb.core import *

marpodb = MarpoDB("/marpodb4")
gene = marpodb.session.query(Gene).first()

print gene.record
