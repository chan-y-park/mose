# TODO: Consider using ConfigParser
import sympy as sym
from math import pi




# #-----------	PURE SU(2)	 ------------

u = sym.Symbol('u')
g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

### Giving by hand the charges at branch points, 
### must update with algorithm that determines actual charge at branch-point
FIXED_CHARGES = [ [1, 0], [-1, 2]] 

### the angle at which branch-cuts should be oriented
THETA_CUTS = 3.14159265359 / 2  

KS_FILTRATION_DEGREE = 4

# how far way from the singularity the locus of the branch cut extends
BRANCH_CUT_CUTOFF = 10.0         

DSZ_MATRIX = [[0, 1], [-1, 0]]

#theta =   0.0 #2.35619449019   ### the phase of K-wall evolution

# the NINT_RANGE are for numerical integration: 
# [initial time, final time, number of steps]
PRIMARY_NINT_RANGE = [0.0, 0.2, 50]
NINT_RANGE = [0.0, 8.0, 400]

THETA_RANGE = [0,pi,100]

# Number of iterations to construct additional K-walls
N_ITERATIONS = 2
#
#-------------------------------------------


# #-----------   SU(2) N_f=1   -------------
# #-------- doesn't work very well ---------

# u = sym.Symbol('u')
# Lambda = 2.5 #- 0.5 * 1j
# m = 0.0#0.1 

# g2 = - Lambda**3 * m + 4 * u**2 / 3
# g3 = Lambda**6 / 16 - m * u * Lambda**3 / 3 + 8 * u**3 / 27

# ### Giving by hand the charges at branch points, 
# ### must update with algorithm that determines actual charge at branch-point
# FIXED_CHARGES = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] 

# ### the angle at which branch-cuts should be oriented
# THETA_CUTS = 3.14159265359 / 2  

# KS_FILTRATION_DEGREE = 4

# # how far way from the singularity the locus of the branch cut extends
# BRANCH_CUT_CUTOFF = 10.0         

# DSZ_MATRIX = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]

# theta =   1.1 ### the phase of K-wall evolution

# # the NINT_RANGE are for numerical integration: 
# # [initial time, final time, number of steps]
# PRIMARY_NINT_RANGE = [0.0, 0.3, 50]
# NINT_RANGE = [0.0, 3.0, 50]

# THETA_RANGE = [0,pi,30]

# N_ITERATIONS = 0


#-------------------------------------------



#-----------   AN INVENTED FIBRATION   -------------

# u = sym.Symbol('u')
# Lambda = 1.0 #- 0.5 * 1j
# m = 0.0#0.1 

# g2 = ((16 * m**2) / 3 - 12 * u / 3) / Lambda**2
# g3 = ((8 * m * u) / Lambda**3 - (64 * m**3) / (27 * Lambda**3) - 4)

# ### Giving by hand the charges at branch points, 
# ### must update with algorithm that determines actual charge at branch-point
# FIXED_CHARGES = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] 

# ### the angle at which branch-cuts should be oriented
# THETA_CUTS = (3.14159265359/2)  

# KS_FILTRATION_DEGREE = 4

# # how far way from the singularity the locus of the branch cut extends
# BRANCH_CUT_CUTOFF = 10.0         

# DSZ_MATRIX = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]

# theta =   12.0/50 * pi ### the phase of K-wall evolution

# # the NINT_RANGE are for numerical integration: 
# # [initial time, final time, number of steps]
# PRIMARY_NINT_RANGE = [0.0, 0.5, 50]
# NINT_RANGE = [0.0, 10.0, 250]

# THETA_RANGE = [0,pi,100]

# N_ITERATIONS = 2


#-------------------------------------------


#------------	AD N=4 THEORY	 ------------
# #-------- doesn't work very well ---------

### Or, really any other theory with a curve  y^2 = x^4 + ...

# u = sym.Symbol('u')
# z = sym.Symbol('z')

# Lambda = 1
# m = 0	### you will get in trouble if m is not zero, 
		### the problem is getting quartic solutions
# P = z**4 + 4 * Lambda**2 * z**2 + 2 * m * z + u
# e1, e2, e3, e4 = solve(P, z)#[lambdify(u, x) for x in solve(P, z)]

# q = simplify(((e2 - e3) * (e4 - e1)) / ((e2 - e1) * (e4 - e3)))

# # print q

# g2 = simplify((4.0 / 3) * (q**2 - q + 1))
# g3 = simplify((4.0 / 27) * (2 * q**3 - 3 * q**2 - 3 * q + 2)) 

# # print g2
# # print g3

# ### Giving by hand the charges at branch points, 
# ### must update with algorithm that determines actual charge at branch-point
# FIXED_CHARGES = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] 

# ### the angle at which branch-cuts should be oriented
# THETA_CUTS = 3.14159265359 / 2  

# KS_FILTRATION_DEGREE = 4

# # how far way from the singularity the locus of the branch cut extends
# BRANCH_CUT_CUTOFF = 10.0         

# DSZ_MATRIX = [[0, 1, -1], [-1, 0, 0], [1, 0, 0]]

# theta =   0.0 #2.35619449019   ### the phase of K-wall evolution

# # the NINT_RANGE are for numerical integration: 
# # [initial time, final time, number of steps]
# PRIMARY_NINT_RANGE = [0.0, 1.0, 101]
# NINT_RANGE = [0.0, 2.0, 400]

# THETA_RANGE = [0,pi,100]

# N_ITERATIONS = 2


#-------------------------------------------


# The following parameter controls what value of abs(det(pf_matrix)) wilL
# raise an exception to determine that a singularity ran too close to a 
# singularity, and should be dropped.
TRAJECTORY_SINGULARITY_THRESHOLD = 10**6 

# The range on the moduli space to search for intersections between 
# K-walls. By default this is also the plot range.
INTERSECTION_SEARCH_RANGE = [[-10, 10], [-10, 10]]
# Size of a bin for the coarse-graining of the intersection module
INTERSECTION_SEARCH_BIN_SIZE = .2 

PF_ODEINT_MXSTEP = 5000000

### Options for saving to files
WRITE_TO_FILE = False
PICKLE_PROTOCOL = 0
### Available options: 0, 1, 2 (up to python 2.7.8)

