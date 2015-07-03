import numpy as np

def charge_monodromy(charge, branch_point, direction):
	"""
	Computes the change in the local charge of a k-wall 
	across a branch cut.
	Given the initial charge gamma, it returns the new
	charge gamma . M
	Monodromy matrices act FROM THE RIGHT in our conventions,
	we follow Seiberg-Wittten.
	"""
	gamma = np.array(charge)
	if direction == 'ccw':
		monodromy_matrix = np.array(branch_point.monodromy_matrix)
		return map(int, map(round, np.dot(charge, monodromy_matrix)))
	if direction == 'cw':
		monodromy_matrix = \
						np.linalg.inv(np.array(branch_point.monodromy_matrix))
		return map(int, map(round, np.dot(charge, monodromy_matrix)))


def flavor_charge_monodromy(gauge_charge, flavor_charge, branch_point, 
															direction, dsz):
	
	"""
	The monodromy matrix of a branch point with 
	gauge charge [p, q] and flavor charge [0,...,1,...,0]
	is given in blocks:
		
			M 	|	A 
		   -----------
			0	|	1
	
	where M is th gauge charge monodromy, A is a 2 x N_f matrix
	witht eh following entries: 
	
			0	0	... 	0	-q 	0	...
			0	0	... 	0	p 	0	...
	
	in such a way that the (p, q, 0, ..., 1,... 0) charge of a branch
	point is a left eigenevctor for the monodromy.
	So the monodromy of a certain vector (m, n, 0, ... , 1, ..., 0) 
	will be (m, n) -> (m, n).M while 
	(0, ..., 1, ..., 0) --> (0, ..., 1, ..., <(m,n), (p,q)> ,..., 0)
	
	Recalling that by default we assign a certain flavor charge to a 
	certain branch point, which is of the form:
	[0,...,1,...,0]
	we can express this in the following manner.
	"""

	bp_gauge_charge = np.array(branch_point.gauge_charge)
	bp_flavor_charge = np.array(branch_point.flavor_charge)

	k_wall_gauge_charge = np.array(gauge_charge)
	k_wall_flavor_charge = np.array(flavor_charge)

	dsz_matrix = np.array(dsz)
	pairing = k_wall_gauge_charge.dot(dsz_matrix.dot(bp_gauge_charge))

	return map(int, map(round, list(k_wall_flavor_charge \
												+ pairing * bp_flavor_charge))) 




