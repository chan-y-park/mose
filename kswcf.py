def progeny_2(data):
	#### The formatting of "data" should be as follows:
	#### [ [gamma_1 , omega_1]  ,  [gamma_2 , omega_2] ] 
	#### phase ordered from right to left
	gamma_1 = data[0][0]
	gamma_2 = data[1][0]
	omega_1 = data[0][1]
	omega_2 = data[1][1]

	# Insert actual KSWCF algorithm ....

	return [ \
	# [[1, 0], 1] ,\
	 [[2, 1], 1] , [[1, 1], -2], [[1,2], 1] \
	#, [[0,1], 1] \
	]
