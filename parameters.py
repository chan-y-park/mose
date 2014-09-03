import sympy as sym
from sympy import *
from sympy import mpmath as mp


u = sym.Symbol('u')
g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

theta = N(pi/2) # N(sym.pi) / 2 

# the options are for numerical integration: initial time, final time, number of steps
primary_options = [0.0, 1.0, 100]
options = [0.0, 3, 100]


# how far way from the singularity the locus of the branch cut extends
branch_cut_cutoff = 10.0         



#### this parameter is related only to the temporary intersection algorithm:
#### this is the threshold above which a distance between two points is not regarded as intersection
intersection_range = 0.1 

dsz_matrix = [[0, 1], [-1, 0]]

new_kwalls = []
kwalls = []
intersections = []

n_iter = 2