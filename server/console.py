from partsdb.partsdb import PartsDB
from tables import *
import os
from partsdb.tools.Exporters import GenBankExporter

if "MARPODB_DB_NAME" in os.environ:
	dbName = '/' + os.environ["MARPODB_DB_NAME"]
else:
	dbName = "marpodb.io/marpodb3?user=common"

marpodb = PartsDB('postgresql://' + dbName, Base = Base)
session = marpodb.Session()
exporter = GenBankExporter(None)
