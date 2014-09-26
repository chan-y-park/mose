"""
Main module. Run 'python -m mose' at the directory that contains the 
directory 'mose'.
"""
import logging
import sys, getopt

from parameters import g2, g3, THETA_CUTS, THETA_RANGE, N_ITERATIONS
from elliptic_fibration import EllipticFibration
from construction import construct_k_wall_networks 
from plotting import ms_plot

try:
    opts, args = getopt.getopt(argv[1:], 'm:', 
                                ['debug', 'info', 'warning']) 
except getopt.GetOptError:
    print 'Unknown options.'

    for opt, arg in opts:
        if (opt == '-m'and arg == 'debug') or opt == '--debug':
            logging_level = logging.DEBUG
            logging_format='%(module)s@%(lineno)d: %(message)s'
        elif (opt == '-m'and arg == 'info') or opt == '--info':
            logging_level = logging.INFO
            logging_format='%(message)s'
        elif (opt == '-m'and arg == 'warning') or opt == '--warning':
            logging_level = logging.WARNING
            logging_format='%(message)s'
        else:
            logging_level = logging.WARNING
            logging_format='%(message)s'

logging.basicConfig(level=logging_level, format=logging_format)

fibration = EllipticFibration(g2, g3, THETA_CUTS)

k_wall_networks = construct_k_wall_networks(THETA_RANGE, N_ITERATIONS)

ms_walls = build_ms_walls(k_wall_networks)

ms_plot(ms_walls)
