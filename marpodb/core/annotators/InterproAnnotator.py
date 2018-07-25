import sys
import xml.etree.ElementTree as ET

from partsdb.tools.Annotators import Annotator

class InterproAnnotator(Annotator):

	def addAnnotation(self, db, dbid, origin, description = '', start = 0, end = 0, refID = '', eVal = 0.0):
		if description == '' or refID == '':
			print "Adding empty {0}({1}) from {2} at {3}:{4} with {5}".format(refID, description, origin, start, end, eVal)
		cds = db.session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == dbid).first()
		if cds:
			match = self.cls(target = cds, description = description, origin = origin, start = start, end = end, refID = refID, eVal= eVal)
			db.session.add(match)
			db.commit()
		else:
			raise Exception("Could not find {0} in the database".format(dbid))

	def annotate(self, db, xmlFile):
		tree = ET.parse(xmlFile)
		matches = tree.getroot()

		ns = {"interpro": "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"}

		for i,protein in enumerate(matches):
			dbid = protein.find('interpro:xref', ns).get('id')

			for match in protein.find('interpro:matches', ns):

				# Getting infor from match (eValue)
				eVal = float( match.get('evalue') ) if match.get('evalue') else 0.0

				# Getting info from signature (origin, description, refID)
				signature = match.find('interpro:signature', ns)

				if signature is not None:
					signatureLibrary = signature.find('interpro:signature-library-release', ns)
					origin = "{0}v.{1}".format(signatureLibrary.get('library'), signatureLibrary.get('version'))
					
					description = signature.get('name') if signature.get('name') else ''
					description = signature.get('desc') if signature.get('desc') else description

					refID = signature.get('ac') if signature.get('ac') else ''

					entry = signature.find('interpro:entry', ns)
					if entry is not None:
						description = entry.get('desc') if entry.get('desc') else description
						refID = entry.get('ac')

						for goterm in entry.findall('interpro:go-xref', ns):
							goDescription = "{0}:{1}".format(goterm.get('category'), goterm.get('name'))
							goID = goterm.get('id')
							self.addAnnotation(db, dbid, "GO", description = goDescription, refID = goID)

				else:
					raise Exception("Could not find signature for {0}".format(dbid))

				# Getting info from location (star, end)
				locations = match.find('interpro:locations', ns)

				if locations is not None:
					location = locations[0]

				if location is not None:
					start = location.get('start')
					end   = location.get('end')
				else:
					raise Exception("Could not find location for {0}".format(dbid))
				
				self.addAnnotation(db, dbid, origin, description, start, end, refID, eVal)

if __name__ == "__main__":
	annotator = InterproAnnotator(None)
	annotator.annotate(None, sys.argv[1])