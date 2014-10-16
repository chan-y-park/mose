import numpy as np

def charge_monodromy(charge, branch_cut, direction):
	gamma = np.array(charge)
	if direction == 'cw':
		monodromy_matrix = np.array(branch_cut.monodromy_matrix)
		return np.dot(charge, monodromy_matrix)
	if direction == 'ccw':
		monodromy_matrix = np.linalg.inv(np.array(branch_cut.monodromy_matrix))
		return dot(charge, monodromy_matrix)