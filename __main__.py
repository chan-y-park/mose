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
    THETA_RANGE
)
from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from k_wall import KWall
from marginal_stability_wall import build_ms_walls
from plotting import plot_k_wall_network, plot_ms_walls
from save_to_file import f_save, f_recover, save_k_wall_network_plot, \
                            save_phase_scan
from misc import formatted_date_time

# Default logging
logging_level = logging.WARNING
logging_format='%(message)s'

generate_single_network = False
generate_multiple_networks = False
write_to_file = False
show_graphics = False

try:
    opts, args = getopt.getopt(sys.argv[1:], 'l:s:fwg', 
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
        Will also produce saved pictures.
        
    -g:
        Show plot of the computation.
        If working at a fixed phase (option -s), it will 
        show a plot of the K-wall network.
        If working with multiple phases (option -f), it
        will show the plot of marginal stability walls.
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
            write_to_file = True
        
        if opt == '-g':
            # save data to external file
            show_graphics = True

except getopt.GetoptError:
    print 'Unknown options.'

logging.basicConfig(level=logging_level, format=logging_format)

fibration = EllipticFibration(g2, g3, FIXED_CHARGES,
                                THETA_CUTS, BRANCH_CUT_CUTOFF)

KWall.count = 0
KWallNetwork.count = 0

### Do not move, must be here for ocnsistency of file naming.
if write_to_file:
    date_time = formatted_date_time()

if generate_single_network == True:
    kwn = KWallNetwork(phase, fibration, INTERSECTION_SEARCH_RANGE,
                        INTERSECTION_SEARCH_BIN_SIZE)
    kwn.grow(PRIMARY_NINT_RANGE, NINT_RANGE, N_ITERATIONS)
    if write_to_file:
        ### save picture
        file_name = 'single_network_' + date_time + '.png'
        save_k_wall_network_plot(kwn, file_name)
        ### save kwn data
        file_name = 'single_network_' + date_time + '.mose'
        saved = f_save(kwn, file_name)
        print saved
    if show_graphics:
        plot_k_wall_network(kwn) 

elif generate_multiple_networks == True:
    k_wall_networks = construct_k_wall_networks(
        fibration, INTERSECTION_SEARCH_RANGE, INTERSECTION_SEARCH_BIN_SIZE,
        PRIMARY_NINT_RANGE, NINT_RANGE, N_ITERATIONS,
        THETA_RANGE
    )
    ms_walls = build_ms_walls(k_wall_networks)
     
    if write_to_file:
        ### save pictures
        file_name_part = 'phase_scan_' + date_time
        save_phase_scan(k_wall_networks, ms_walls, file_name_part, \
                                                INTERSECTION_SEARCH_RANGE)
        ### save all data
        file_name = 'phase_scan_' + date_time + '.mose'
        saved = f_save([k_wall_networks, ms_walls], file_name)
        print saved
    if show_graphics:
        plot_ms_walls(ms_walls, INTERSECTION_SEARCH_RANGE)
