import logging

from structure import build_first_generation, iterate, prepare_branch_locus 
from trajectory_animation import animate_trajectories
from parameters import *

logging.basicConfig(level=logging.DEBUG, format=LOGGING_FORMAT)
#logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.WARNING)

new_kwalls = []
kwalls = []
intersections = []

logging.info('Preparing the branch locus:')

branch_locus = prepare_branch_locus(g2, g3, theta)
branch_points = branch_locus[0]
branch_cuts = branch_locus[1]

logging.info('Branch point loci and charges: %s ', 
                [[branch_points[i].locus, branch_points[i].charge] 
                    for i in range(len(branch_points))]
)

logging.info('Branch cuts loci and charges: %s ', 
                [[branch_cuts[i].locus, branch_cuts[i].charge] 
                    for i in range(len(branch_cuts))]
)

logging.info('Preparing the primary kwalls:')

new_kwalls = build_first_generation(branch_points, theta, g2, g3, 
                                    primary_options, options)


logging.info('Iterating %s times:', n_iter)

kwalls, new_kwalls, intersections = iterate(n_iter, kwalls, new_kwalls, 
                                            intersections)

animate_trajectories(kwalls+new_kwalls, 1)

# phase_scan(theta_range)
