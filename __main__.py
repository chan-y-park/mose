"""
Main module. Run 'python -m mose' at the directory that contains the 
directory 'mose'.
"""

import logging

from mose import (branch_points, theta, g2, g3, n_iter)
from structure import build_first_generation, iterate, prepare_branch_locus 
from trajectory_animation import animate_trajectories

new_kwalls = []
kwalls = []
intersections = []

logging.info('Preparing the primary kwalls:')

new_kwalls = build_first_generation(branch_points, theta, g2, g3)

logging.info('Iterating %s times:', n_iter)

kwalls, new_kwalls, intersections = iterate(n_iter, kwalls, new_kwalls, 
                                            intersections)

animate_trajectories(kwalls+new_kwalls, 1)

# phase_scan(theta_range)
