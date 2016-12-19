def compRev(seq, sign):
	if sign =='-':
		sub = dict( zip( ['A', 'T', 'G', 'C', 'N'], ['T', 'A', 'C', 'G', 'N'] ) )
		return ''.join( sub[nucl] for nucl in seq.upper() )[::-1]
	else:
		return seq