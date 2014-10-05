"""
Main module. Run 'python -m mose' at the directory that contains the 
directory 'mose'.
"""
import logging
import sys, getopt

from config import (
    g2, g3, 
    FIXED_CHARGES, THETA_CUTS, BRANCH_CUT_CUTOFF, N_ITERATIONS, 
    PRIMARY_NINT_RANGE, NINT_RANGE,
    INTERSECTION_SEARCH_RANGE, INTERSECTION_SEARCH_BIN_SIZE,
    THETA_RANGE, WRITE_TO_FILE
)
from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from marginal_stability_wall import build_ms_walls
from plotting import plot_k_wall_network, ms_plot
from misc import formatted_date_time
from save_to_file import f_save, f_recover

# Default logging
logging_level = logging.WARNING
logging_format='%(message)s'

generate_single_network = False
generate_multiple_networks = False

try:
    opts, args = getopt.getopt(sys.argv[1:], 'l:s:f:w', 
                                ['logging_level=']) 

    if len(opts) == 0:
        print("""usage: python -m mose [OPTION]

    -l LEVEL, --logging_level=LEVEL:
        set logging level to LEVEL.

    -s THETA:
        produce a single K-wall network at the phase of THETA
        and plot the K-wall network.

    -f:
        produce K-wall networks and plot walls of marginal
        stability.

    -w:
        save data on a date-named file for offline analysis,
        the name will include 'phase_scan' or 'single_network'
        according to what data is stored.
        """
        )

    for opt, arg in opts:
        if (opt == '-l' or opt == '--logging_level'):
            if arg == 'debug':
                logging_level = logging.DEBUG
                logging_format='%(module)s@%(lineno)d: %(message)s'
            elif arg == 'info':
                logging_level = logging.INFO
            elif arg == 'warning':
                logging_level = logging.WARNING

        if opt == '-s':
            # Generate a single K-wall network at a phase
            phase = float(arg)
            generate_single_network = True 
        elif opt == '-f':
            # Generate K-wall networks at various phases
            generate_multiple_networks = True

        if opt == '-w':
            # save data to external file
            WRITE_TO_FILE = True

except getopt.GetoptError:
    print 'Unknown options.'

logging.basicConfig(level=logging_level, format=logging_format)

fibration = EllipticFibration(g2, g3, FIXED_CHARGES,
                                THETA_CUTS, BRANCH_CUT_CUTOFF)

if generate_single_network == True:
    kwn = KWallNetwork(phase, fibration, INTERSECTION_SEARCH_RANGE,
                        INTERSECTION_SEARCH_BIN_SIZE)
    kwn.grow(PRIMARY_NINT_RANGE, NINT_RANGE, N_ITERATIONS)
    if WRITE_TO_FILE:
        file_name = 'single_network_' + formatted_date_time() + '.mose'
        saved = f_save(kwn, file_name)
        print saved
    plot_k_wall_network(kwn) 

elif generate_multiple_networks == True:
    k_wall_networks = construct_k_wall_networks(
        fibration, INTERSECTION_SEARCH_RANGE, INTERSECTION_SEARCH_BIN_SIZE,
        PRIMARY_NINT_RANGE, NINT_RANGE, N_ITERATIONS,
        THETA_RANGE
    )
    ms_walls = build_ms_walls(k_wall_networks)
    ms_plot(ms_walls, INTERSECTION_SEARCH_RANGE)
