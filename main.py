from numerics import *
from structure import *
from trajectory_animation import animate_trajectories
from parameters import *





new_kwalls = []
kwalls = []
intersections = []

if verb: print "Preparing the branch locus:"
branch_locus = prepare_branch_locus(g2, g3, theta)
branch_points = branch_locus[0]
branch_cuts = branch_locus[1]
if verb: print "Branch point loci and charges: %s " % [ [branch_points[i].locus, branch_points[i].charge] for i in range(len(branch_points)) ]
if verb: print "Branch cuts loci and charges: %s " % [ [branch_cuts[i].locus, branch_cuts[i].charge] for i in range(len(branch_cuts)) ]

if verb: print "\nPreparing the primary kwalls:"
new_kwalls = build_first_generation(branch_points, theta, g2, g3, primary_options, options)


if verb: print "\nIterating %s times:" % n_iter
kwalls, new_kwalls, intersections = iterate(n_iter, kwalls, new_kwalls, intersections)

animate_trajectories(kwalls+new_kwalls,1)


# phase_scan(theta_range)