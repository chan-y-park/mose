"""
Main module. Run 'python -m mose' at the directory that contains the
directory 'mose'.
"""
import logging
import os
import sys
import getopt

from mose_config_parser import MoseConfigParser
from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from k_wall import KWall
from marginal_stability_wall import build_ms_walls
from plotting import plot_k_wall_network, plot_ms_walls
from save_to_file import (f_save, f_recover, save_k_wall_network_plot,
                          save_phase_scan, prepare_folder)
from misc import formatted_date_time

config = MoseConfigParser()
config_file = ''

# Default logging
logging_level = logging.WARNING
logging_format = '%(message)s'

generate_single_network = False
generate_multiple_networks = False
write_to_file = False
show_graphics = False

try:
    opts, args = getopt.getopt(sys.argv[1:], 'c:l:s:fwg', ['logging_level=']) 

    if len(opts) == 0:
        print("""usage: python -m mose [OPTION]

    -c CFG_FILE_NAME:
        read CFG_FILE_NAME to set up the configuration,
        including fibration, ODE solver parameters, etc.

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
        """)

    for opt, arg in opts:
        if (opt == '-c' and len(arg) > 0):
            config_file = arg
        if (opt == '-l' or opt == '--logging_level'):
            if arg == 'debug':
                logging_level = logging.DEBUG
                logging_format = '%(module)s@%(lineno)d: %(message)s'
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



if generate_single_network or generate_multiple_networks:
    logging.basicConfig(level=logging_level, format=logging_format)

    main_file_dir, main_file_name = os.path.split(__file__)
    if len(config_file) == 0:
        # default configuration file
        config_file = 'fibration_invented.ini'
        logging.warning('No .ini file specified --- load %s instead.', config_file)
    config.read(os.path.join(main_file_dir, config_file))
    logging.debug('Configuration sections: %s', config.sections())
    logging.debug('g2 = %s', config.get('fibration', 'g2'))
    logging.debug('g3 = %s', config.get('fibration', 'g3'))

    fibration = EllipticFibration(
        config.get('fibration', 'g2'),
        config.get('fibration', 'g3'),
        config.get('charge', 'fixed_charges'),
        config.get('charge', 'dsz_matrix'),
        config.get('branch cut', 'theta'),
        config.get('branch cut', 'cutoff')
    )

    KWall.count = 0
    KWallNetwork.count = 0

# DO NOT MOVE: must be here for consistency of file naming.
if write_to_file:
    if generate_single_network:
        label = 'single_network'
    elif generate_multiple_networks:
        label = 'phase_scan'
    current_dir = os.getcwd()
    date_time = formatted_date_time()
    plots_dir = prepare_folder(label + '_' + date_time)

if generate_single_network is True:
    kwn = KWallNetwork(
        phase, fibration,
        config.get('intersection search', 'range'),
        config.get('intersection search', 'bin_size')
    )
    kwn.grow(
        config.get('ODE', 'primary_k_wall_odeint_range'),
        config.get('ODE', 'odeint_range'),
        config.get('ODE', 'trajectory_singularity_threshold'),
        config.get('ODE', 'pf_odeint_mxstep'),
        config.get('K-wall network', 'n_iterations'),
        config.get('KSWCF', 'filtration_degree')
    )
    if write_to_file:
        # save picture
        file_path = os.path.join(plots_dir, 
                        'single_network_' + date_time + '.png')
        save_k_wall_network_plot(kwn, file_path)
        # save kwn data
        file_path = os.path.join(current_dir, 'results', 
                        'single_network_' + date_time + '.mose')
        saved = f_save(kwn, file_path, config.get('file IO', 
                        'pickle_protocol'))
        print saved
    if show_graphics:
        plot_k_wall_network(kwn) 

elif generate_multiple_networks is True:
    k_wall_networks = construct_k_wall_networks(
        fibration,
        config.get('intersection search', 'range'),
        config.get('intersection search', 'bin_size'),
        config.get('ODE', 'primary_k_wall_odeint_range'),
        config.get('ODE', 'odeint_range'),
        config.get('ODE', 'trajectory_singularity_threshold'),
        config.get('ODE', 'pf_odeint_mxstep'),
        config.get('K-wall network', 'n_iterations'),
        config.get('KSWCF', 'filtration_degree'),
        config.get('MS wall', 'theta_range')
    )
    ms_walls = build_ms_walls(k_wall_networks)
     
    if write_to_file:
        # save pictures
        file_path_part = os.path.join(plots_dir, 'phase_scan_' + date_time)
        save_phase_scan(k_wall_networks, ms_walls, file_path_part,
                        config.get('intersection search', 'range'))
        # save all data
        file_path = os.path.join(current_dir, 'results', 
                        'phase_scan_' + date_time + '.mose')
        saved = f_save([k_wall_networks, ms_walls], file_path,
                       config.get('file IO', 'pickle_protocol'))
        print saved
    if show_graphics:
        plot_ms_walls(ms_walls, config.get('intersection search', 'range'))
