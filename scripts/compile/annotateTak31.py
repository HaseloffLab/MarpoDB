from marpodb.core import *
from sqlalchemy import and_

# Use Info_names list from marchantia.info to annotate CDSs
def annotateInfoNames(client, infoNames = "../../data/marpodb4/Info_names"):
	
	session = client.session

	infoFile = open(infoNames)
	for line in infoFile:
		tabs = line.rstrip().split()

		print line

		name = tabs[0]
		alias = tabs[1]

		dbxref = session.query(DbxRef).filter(DbxRef.refID == alias).first()

		if dbxref:
			gene = dbxref.target.gene[0]
			
			gene.name = name
			dbref = session.query(DbxRef).filter(and_(DbxRef.origin == "MarpolBase", DbxRef.targetID == gene.cds.id ) ).first()
			if dbref:
				if dbref.description != name:
					print "Changing {0} name from {1} to {2}".format(gene.cds.dbid, dbref.description, name)
					dbref.description = name
			else:
				client.addPart("dbxref", description = name, refID = '', origin = "MarpolBase", target = gene.cds)		
		else:
			print "Didn't find ", alias

	client.commit()   		

	session.close()

if __name__ == "__main__":
	marpodb = MarpoDB("/marpodb4")
	# marpodb.annotate("interprohit", xmlFile = "../../data/marpodb4/marpodb4_cam_tak_2018_06_27.fa.xml")
	annotateInfoNames(marpodb)
	# marpodb.annotate("blastphit", blastFileName = "../../data/marpodb4/marpodb4_cam_tak_2018_06_27_vs_uniprot_2018_07_02.outfmt6", uniprotInfoFileName = "../../data/misc/uniprot/uniprot_viridiplantae_trans_prot_evidence_02_07_2018.tab")
	