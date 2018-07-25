from marpodb.core import *

marpodb = MarpoDB("/marpodb4")

session = marpodb.Session()

cdss = session.query(CDS).all()
marpodb.export(cdss, "../../output/state/CDS.fa", pep = True)