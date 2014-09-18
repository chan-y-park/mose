import logging
import sympy as sym

from math import pi
from structure import build_first_generation, iterate, prepare_branch_locus 
from trajectory_animation import animate_trajectories

#---------------------------------------------------------------------------
# Program options
#---------------------------------------------------------------------------

LOGGING_FORMAT = '%(levelname)s:%(module)s@(%(lineno)d): %(message)s'
LOGGING_LEVEL = logging.DEBUG
#LOGGING_LEVEL = logging.INFO
#LOGGING_LEVEL = logging.WARNING

logging.basicConfig(level=LOGGING_LEVEL, format=LOGGING_FORMAT)

#---------------------------------------------------------------------------
# Default parameters
#---------------------------------------------------------------------------

u = sym.Symbol('u')
g2 = 1 + 4 * (u ** 2) / 3
g3 = u / 3 + 8 * (u ** 3) / 27

theta = pi 

# the options are for numerical integration: initial time, final time, number 
# of steps.
primary_options = [0.0, 1.0, 100]
options = [0.0, 1.0, 100]

# how far way from the singularity the locus of the branch cut extends
branch_cut_cutoff = 10.0         


NUMBER_OF_ITERATIONS = 2

THETA_I = 0
THETA_F = pi
DELTA_THETA = 3

THETA_RANGE = [THETA_I, THETA_F, DELTA_THETA] 

KS_FILTRATION_DEGREE = 6

#NOTE: Add a routine that accepts sys argv and set parameters accordingly.
if True:
    n_iter = NUMBER_OF_ITERATIONS

#---------------------------------------------------------------------------
# Objects
#---------------------------------------------------------------------------

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
