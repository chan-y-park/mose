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
	if direction == 'cw':
		monodromy_matrix = np.array(branch_point.monodromy_matrix)
		return map(int, map(round, np.dot(monodromy_matrix, charge)))
	if direction == 'ccw':
		monodromy_matrix = \
						np.linalg.inv(np.array(branch_point.monodromy_matrix))
		return map(int, map(round, np.dot(monodromy_matrix, charge)))