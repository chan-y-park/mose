import numpy as np

def charge_monodromy(charge, branch_point, direction):
	gamma = np.array(charge)
	if direction == 'cw':
		monodromy_matrix = np.array(branch_point.monodromy_matrix)
		return map(int, map(round, np.dot(charge, monodromy_matrix)))
	if direction == 'ccw':
		monodromy_matrix = \
						np.linalg.inv(np.array(branch_point.monodromy_matrix))
		return map(int, map(round, np.dot(charge, monodromy_matrix)))