# Add dataset to marpodb3 records
import sys

from marpodb.core import *

marpodb = MarpoDB("/marpodb4")
marpodb.setup(prefix="mpdb")
Base.metadata.create_all(marpodb.engine)

datasetName = "cam"
datasetVersion = "1.0"

dataset = marpodb.session.query(Dataset).filter(Dataset.name == datasetName, Dataset.version == datasetVersion).first()
if not dataset:
	dataset = marpodb.addPart("dataset", name = datasetName, version = datasetVersion)

genes = marpodb.session.query(Gene).all()

for gene in genes:
	gene.dataset = dataset

marpodb.commit()
