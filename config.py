#TODO: Consider using ConfigParser
import sympy as sym
from math import pi

u = sym.Symbol('u')
g2 = 1 + 4*(u ** 2)/3
g3 = u/3 + 8*(u ** 3)/27

# The angle at which branch-cuts should be oriented.
THETA_CUTS = (3.14159265359/2)  

# The filtration of KSWCF
KS_FILTRATION_DEGREE = 4

# How far way from the singularity the locus of the branch cut extends.
BRANCH_CUT_CUTOFF = 10.0         

DSZ_MATRIX = [[0, 1], [-1, 0]]


# The following parameter controls what value of abs(det(pf_matrix)) wilL
# raise an exception to determine that a singularity ran too close to a 
# singularity, and should be dropped.
TRAJECTORY_SINGULARITY_THRESHOLD = 10**6 

# The range on the moduli space to search for intersections between 
# K-walls. By default this is also the plot range.
INTERSECTION_SEARCH_RANGE = [[-10, 10], [-10, 10]]
# Size of a bin for the coarse-graining of the intersection module
INTERSECTION_SEARCH_BIN_SIZE = .2 

# Options for numerical integration: 
# [initial time, final time, number of steps]
PRIMARY_NINT_RANGE = [0.0, 1.0, 101]
NINT_RANGE = [0.0, 2.0, 400]

# Range of the phase to scan:
# [initial phase, final phase, number of steps]
THETA_RANGE = [0, pi, 8]

# Number of iterations to construct additional K-walls
N_ITERATIONS = 2

# TODO: Must update with actual charge at branch-point
FIXED_CHARGES = [[1, 0], [-1, 2]] 

PF_ODEINT_MXSTEP = 5000000
