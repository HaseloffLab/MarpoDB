# MarpoDB 4
This is a stripped down version of the master branch, which is enough to run a local version of the MarpoDB, both backend and frontend.

## Requirements
* PostgreSQL >= 9.5.2
* Python 2.7
* virtualenv

## Setting up

We recommend to use virtual environemt. You can create on by
``` bash
virtualenv -p python2.7 marpodb
```
and activate it with
``` bash
source marpodb/bin/activate
```

### Automatic

1. Clone the marpodb4 branch to you virtual environment folder:
``` bash
git clone https://github.com/HaseloffLab/MarpoDB.git -b marpodb4 marpodb4
```
2. Install the package using pip:
```bash
pip install .
```
3. This should install the MarpoDB package toogether with an installation script, which will create a database for you, download the dump from the server and restore it locally. To do so, simply run
```bash
installmarpodb -dbname marpodb
```

### Manual
```
2. Create a new PostgreSQL database:
``` bash
createdb marpodb
```
3. Download the dump of the database from marpodb.io/data, and restore it into the newly-created database:
``` bash
pg_restore -d marpodb marpodb4.dump
```
this might take a while.

4. Clone the marpodb4 branch to you virtual environment folder:
``` bash
git clone https://github.com/HaseloffLab/MarpoDB.git -b marpodb4 marpodb4
```
5. Go to the marpodb4 folder and install the package:
``` bash
pip install .
```
Now you should be ready to access the database using the Python bindings.

## Creating a session

1. The followig code connects to the marpodb database, and creates a SQLAlchemy session object:
``` python
from marpodb.core import *

marpodb = MarpoDB('/marpodb')

session = marpodb.Session()
```
The argument you pass to MarpoDB should be the same as the name of the databse you have created.

2. You can now run queries, for example:
``` python
myGene = session.query(Gene).filter(Gene.dbid == 'mpdb.gene.8238').first()
```
Selects a gene by its dbid.

or:
``` python
myPromoter = session.query(Promoter).filter(Promoter.id == Gene.promoterID)\
.filter(Gene.datasetID == Dataset.id).filter(Dataset.name == 'tak').filter(CDS.id == Gene.cdsID)\
.filter(CDS.id == DbxRef.targetID).filter(DbxRef.refID == 'Mapoly0063s0040').first()
```
You can check out the database schema at marpodb.io/data and learn more about SQLAlchemy.

## Exporting sequences

Once you have selected part or list of parts you export them using a call to the MarpoDB object, e.g.
``` python
marpodb.export( [myPromoter, myGene], "parts.fa", format = "fasta", pep = False)
```
will export the selected parts into "parts.fa" file. The pep argument controls, whether to export CDS as a nucleotide or protein sequence.
