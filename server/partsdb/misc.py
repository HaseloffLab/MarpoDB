
def nonRemoveIDGenerator(cls, session):
	return session.query(cls).count()