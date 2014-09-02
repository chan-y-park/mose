from numerics import *
from structure import *
from trajectory_animation import animate_trajectories

import sympy as sym
from sympy import *
from sympy import mpmath as mp

u = sym.Symbol('u')

g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

theta = N(pi/2) # N(sym.pi) / 2 

# the options are for numerical integration: initial time, final time, number of steps
primary_options = [0.0, 1.0, 100]
options = [0.0, 4, 400]


print "Preparing the branch locus:"
branch_locus = prepare_branch_locus(g2, g3, theta)
branch_points = branch_locus[0]
branch_cuts = branch_locus[1]
print "Branch point loci and charges: %s " % [ [branch_points[i].locus, branch_points[i].charge] for i in range(len(branch_points)) ]
print "Branch cuts loci and charges: %s " % [ [branch_cuts[i].locus, branch_cuts[i].charge] for i in range(len(branch_cuts)) ]

print "\nPreparing the primary kwalls:"
new_kwalls = []
kwalls = []
intersections = []
new_kwalls = build_first_generation(branch_points, theta, g2, g3, primary_options, options)


animate_trajectories(new_kwalls,3)



