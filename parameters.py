#TODO: parameters that are constant should be moved to __init__, parameters 
# that can be changed and thus will be obtained via command line options
# should be moved to __main__ 
import sympy as sym
from sympy import *
from sympy import mpmath as mp
from sympy.solvers import solve
import numpy as np

#TODO: These parameters should be moved to __init__

kwalls = []
new_kwalls = []
intersections = []

# #-----------	PURE SU(2)	 ------------
# u = sym.Symbol('u')
# g2 = 1 + 4 * (u ** 2) / 3
# g3 = u / 3 + 8 * (u ** 3) / 27

# ### Giving by hand the charges at branch points, 
# ### must update with algorithm that determines actual charge at branch-point
# fixed_charges = [ [1, 0], [-1, 2] ] 

# ### the angle at which branch-cuts should be oriented
# theta_cuts = 3.14159265359 / 2  

# ks_filtration_degree = 4

# # how far way from the singularity the locus of the branch cut extends
# branch_cut_cutoff = 10.0         

# dsz_matrix = [[0, 1], [-1, 0]]

#-------------------------------------------


# #-----------   SU(2) N_f=1   -------------
u = sym.Symbol('u')
Lambda = 1
m = 0	### you will get in trouble if m is not zero, the problem is getting quartic solutions
g2 = ((16 * m**2) / 3 - 12 * u / 3) / Lambda**2
g3 = ((8 * m * u) / Lambda**3 - (64 * m**3) / (27 * Lambda**3) - 4)

### Giving by hand the charges at branch points, 
### must update with algorithm that determines actual charge at branch-point
fixed_charges = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] 

### the angle at which branch-cuts should be oriented
theta_cuts = 3.14159265359 / 2  

ks_filtration_degree = 4

# how far way from the singularity the locus of the branch cut extends
branch_cut_cutoff = 10.0         

dsz_matrix = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]

#-------------------------------------------



#------------	AD N=4 THEORY	 ------------
### Or, really any other theory with a curve  y^2 = x^4 + ...

# u = sym.Symbol('u')
# z = sym.Symbol('z')

# Lambda = 1
# m = 0	### you will get in trouble if m is not zero, the problem is getting quartic solutions
# P = z**4 + 4 * Lambda**2 * z**2 + 2 * m * z + u
# e1, e2, e3, e4 = solve(P, z)#[lambdify(u, x) for x in solve(P, z)]

# #----------- Alternative fibration that I invented ----------
# # e1 = u
# # e2 = 2**(-1*u)
# # e3 = 2*(2**(-1*u) + 2**(1*u))
# # e4 = 10*u - 10 + 4**(1*u)
# #-------------------------------------------------------------

# q = simplify(((e2 - e3) * (e4 - e1)) / ((e2 - e1) * (e4 - e3)))

# # print q

# g2 = simplify((4.0 / 3) * (q**2 - q + 1))
# g3 = simplify((4.0 / 27) * (2 * q**3 - 3 * q**2 - 3 * q + 2)) 

# # print g2
# # print g3

# ### Giving by hand the charges at branch points, 
# ### must update with algorithm that determines actual charge at branch-point
# fixed_charges = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ] 

# ### the angle at which branch-cuts should be oriented
# theta_cuts = 3.14159265359 / 2  

# ks_filtration_degree = 4

# # how far way from the singularity the locus of the branch cut extends
# branch_cut_cutoff = 10.0         

# dsz_matrix = [[0, 1, -1], [-1, 0, 0], [1, 0, 0]]

#-------------------------------------------



#-------------------------------------------

### The following parameter controls what value of abs(det(pf_matrix)) will
### raise an exception to determine that a singularity ran too close to a 
### singularity, and should be dropped.
TRAJECTORY_SINGULARITY_THRESHOLD = 10 ** 6 

INTERSECTION_SEARCH_RANGE = [[-10, 10], [-10, 10]]
INTERSECTION_SEARCH_BIN_SIZE = .2 
#TODO: End of __init__ part

#TODO: These parameters should be moved to __main__

#TODO: End of __main__ part

theta =   0.0 #2.35619449019   ### the phase of K-wall evolution

# the options are for numerical integration: 
# [initial time, final time, number of steps]
primary_options = [0.0, 1.0, 101]
options = [0.0, 2.0, 400]

theta_range = [0,np.pi,100]

n_iter = 2

