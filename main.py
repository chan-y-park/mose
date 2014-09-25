import logging

from structure import build_first_generation, iterate, prepare_branch_locus, phase_scan
# from trajectory_animation import animate_trajectories
from parameters import *
from msplot import ms_plot, kwallplot

LOGGING_FORMAT = '%(module)s@%(lineno)d: %(message)s'
#logging.basicConfig(level=logging.DEBUG, format=LOGGING_FORMAT)
logging.basicConfig(level=logging.INFO, format='%(message)s')
#logging.basicConfig(level=logging.WARNING)



logging.info('Preparing the branch locus:')
branch_locus = prepare_branch_locus(g2, g3, theta_cuts)
branch_points = branch_locus[0]
branch_cuts = branch_locus[1]

for i in range(len(branch_points)):
	logging.info('\nBranch point number %s:\n\tlocus %s\n\tcharge %s', 
		branch_points[i].count, branch_points[i].locus, branch_points[i].charge
	)

logging.info('\nPreparing the primary kwalls:')
new_kwalls = build_first_generation(branch_points, theta, g2, g3)


logging.info('\nIterating %s times:', n_iter)
kwalls, new_kwalls, intersections = iterate(n_iter, kwalls, new_kwalls, 
                                          intersections)

# animate_trajectories(kwalls+new_kwalls, 1)

# kwallplot(kwalls+new_kwalls, INTERSECTION_SEARCH_RANGE,
#             branch_points=branch_points) 

# all_data = phase_scan(theta_range)
# all_intersections = all_data[0]
# all_kwalls = all_data[1]
# ms_walls = all_data[2]

# ms_plot(ms_walls)

