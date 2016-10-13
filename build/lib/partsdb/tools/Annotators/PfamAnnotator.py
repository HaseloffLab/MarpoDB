from . import Annotator

class PfamAnnotator(Annotator):
	def annotate(self, db, **kwargs):
		self._annotate(db, kwargs['fileName'])

	def _annotate(self, db, fileName):

		inFile = open(fileName)

		for line in inFile:
			if not line.startswith('#'):
				tabs = line.rstrip().split(' ')
				tabs = filter(None, tabs)

				hit = self.cls()

				hit.name = tabs[0]
				hit.acc = tabs[1]
				cdsID = tabs[3]
				hit.eVal = tabs[6]
				hit.cVal = tabs[11]
				
				locS = tabs[19]
				locE = tabs[20]
				hit.coordinates = str(locS) + ',' + str(locE)

				hit.description = ' '.join(tabs[22:])

				cds = db.session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == cdsID).first()
				if cds:
					hit.target = cds
					# hit.targetID = cds.id
					db.session.add(hit)
				else:
					print 'Failed to locate cds'
		inFile.close()
		db.commit()