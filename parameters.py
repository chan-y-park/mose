import sympy as sym
from sympy import *
from sympy import mpmath as mp
import numpy as np

u = sym.Symbol('u')
g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

theta_cuts = 3.14159265359 / 2  ### the angle at which branch-cuts should be oriented
#theta =   0.942477796077  ### the phase of K-wall evolution
theta =   0  ### the phase of K-wall evolution

# the options are for numerical integration: initial time, final time, number of steps
primary_options = [0.0, 1.0, 100]
options = [0.0, 2.0, 400]


# how far way from the singularity the locus of the branch cut extends
branch_cut_cutoff = 10.0         


#### this parameter is related only to the temporary intersection algorithm:
#### this is the threshold above which a distance between two points is not 
#### regarded as intersection. Drop it once intersection module is ready.
#intersection_range = 0.02 

dsz_matrix = [[0, 1], [-1, 0]]

new_kwalls = []
kwalls = []
intersections = []

n_iter = 2

theta_range = [0,np.pi,10]

verb = True

ks_filtration_degree = 4

INTERSECTION_SEARCH_RANGE = [[-10, 10], [-10, 10]]
INTERSECTION_SEARCH_BIN_SIZE = .2 
#INTERSECTION_SEARCH_BIN_SIZE = .7 
