#TODO: parameters that are constant should be moved to __init__, parameters 
# that can be changed and thus will be obtained via command line options
# should be moved to __main__ 
import sympy as sym
from sympy import *
from sympy import mpmath as mp
import numpy as np

#TODO: These parameters should be moved to __init__

u = sym.Symbol('u')
g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

### the angle at which branch-cuts should be oriented
theta_cuts = 3.14159265359 / 2  

ks_filtration_degree = 4

# how far way from the singularity the locus of the branch cut extends
branch_cut_cutoff = 10.0         

dsz_matrix = [[0, 1], [-1, 0]]


### The following parameter controls what value of abs(det(pf_matrix)) wilL
### raise an exception to determine that a singularity ran too close to a 
### singularity, and should be dropped.
TRAJECTORY_SINGULARITY_THRESHOLD = 10 ** 6 

INTERSECTION_SEARCH_RANGE = [[-10, 10], [-10, 10]]
INTERSECTION_SEARCH_BIN_SIZE = .2 
#TODO: End of __init__ part

#TODO: These parameters should be moved to __main__

#TODO: End of __main__ part

theta =   2.35619449019   ### the phase of K-wall evolution

# the options are for numerical integration: 
# [initial time, final time, number of steps]
primary_options = [0.0, 1.0, 101]
options = [0.0, 2.0, 400]

theta_range = [0,np.pi,100]

n_iter = 2

