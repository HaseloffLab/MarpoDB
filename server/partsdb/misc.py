
def nextIDGenerator(cls, session):
	try:
		ids = [ int(row.id.split('.')[2]) for row in session.query(cls) ]
		return max(ids) + 1
	except:
		return 0
	
