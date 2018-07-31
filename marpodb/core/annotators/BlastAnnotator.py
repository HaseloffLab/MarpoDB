from partsdb.tools.Annotators import Annotator

class BlastAnnotator(Annotator):

	def annotate(self, db, **kwargs):
		if 'blastFileName' in kwargs and 'uniprotInfoFileName' in kwargs:
			self.twoFileAnnotate(db, kwargs['blastFileName'], kwargs['uniprotInfoFileName'])

	def twoFileAnnotate(self, db, blastFileName, uniprotInfoFileName):
		badNames = ('Putative', 'putative', 'Predicted protein', 'Deleted.', 'Merged into', 'Uncharacterized', 'Uncharacterised')
		
		uniprotFile = open(uniprotInfoFileName)

		unidb = {}

		for line in uniprotFile:
			tabs = line.rstrip().split('\t')

			origin = tabs[2].split('(')[0].rstrip()
			proteinName = tabs[1].split('(')[0].rstrip()

			good = True

			for badName in badNames:
				if badName in proteinName:
					good = False
					break
			
			if good:
				unidb[tabs[0]] = [ proteinName, origin ]

		print unidb.keys()
		blastFile = open(blastFileName)
		
		BUFF_SIZE = 1000000
		n = 0

		buff = blastFile.readlines(BUFF_SIZE)
		
		while buff:
			n += 1
			print buff[0] 
			
			# Read hits
			for line in buff:
				tabs = line.rstrip().split('\t')
		
				names = ["cdsID", "uniID", "qlen", "slen", "qstart", "qend", "tstart", "tend", "qcovs", "pident", "evalue"]

				data = dict( zip(names, tabs) )
				
				uniID = data["uniID"].split('|')[1]

				if not uniID in unidb:
					# print "No {0} in unidb".format( uniID )
					continue

				data["uniID"] = uniID
				data["proteinName"] = unidb[uniID][0]
				data["origin"] 		= unidb[uniID][1]
				
				hit = self.cls()
				hit.refID					= data["uniID"]
				hit.coverage				= data["qcovs"]
				hit.qLen	 				= data["qlen"]
				hit.tLen	 				= data["slen"]
				hit.start					= data["qstart"]
				hit.end						= data["qend"]
				hit.tStart					= data["tstart"]
				hit.tEnd					= data["tend"]
				hit.identity				= data["pident"]
				hit.eVal	 				= float(data["evalue"])
				hit.description	 			= data["proteinName"]
				hit.origin	 				= data["origin"]
				hit.eVal 					= float(data["evalue"])
				
				cds = db.session.query(self.cls.__targetclass__).filter(self.cls.__targetclass__.dbid == data['cdsID']).first()

				if cds:
					hit.target = cds
					db.session.add(hit)
				else:
					print "Failed to locate {0}".format(hitName.split('_')[0])
			db.commit()
			buff = blastFile.readlines(BUFF_SIZE)