from partsdb.partsdb import PartsDB
from tables import *
import os
from partsdb.tools.Exporters import GenBankExporter

marpodb = PartsDB('postgresql:///' + os.environ["MARPODB_DB_NAME"], Base = Base)
session = marpodb.Session()
exporter = GenBankExporter(None)
