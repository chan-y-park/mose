import logging
import os
import sys
import getopt
import time

from mose_config_parser import MoseConfigParser
from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from k_wall import KWall
from marginal_stability_wall import build_ms_walls
from plotting import plot_k_wall_network, plot_ms_walls
from save_to_file import (f_save, f_recover, save_k_wall_network_plot,
                          save_phase_scan, prepare_folder)
from misc import formatted_date_time



config_file = ''
config = MoseConfigParser()

def analysis(graphics, save, analysis_type, log, fibration, 
             phase = 0, 
             #theta_range=[0,3.14159,10],
             theta_range=[],
             show_bins = False,
             show_data_points = False,
             show_segments = False,):
    generate_single_network = False
    generate_multiple_networks = False
    write_to_file = False
    show_graphics = False  

    if log == 'debug':
        logging_level = logging.DEBUG
        logging_format = '%(module)s@%(lineno)d: %(message)s'
    elif log == 'info':
        logging_level = logging.INFO
        logging_format = '%(process)d: %(message)s'
    elif log == 'warning':
        logging_level = logging.WARNING
        logging_format = '%(message)s'

    logging.basicConfig(level=logging_level, format=logging_format, 
                        stream=sys.stdout)

    config_file = fibration

    if analysis_type == 'single':
        # Generate a single K-wall network at a phase
        generate_single_network = True
    elif analysis_type == 'full':
        # Generate K-wall networks at various phases
        generate_multiple_networks = True
    if graphics == True:
        # save data to external file
        show_graphics = True
    if save == True:
        # save data to external file
        write_to_file = True

    #-------------------------------------------------------------------------

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    if generate_single_network or generate_multiple_networks:

        main_file_dir, main_file_name = os.path.split(__file__)
        if len(config_file) == 0:
            # default configuration file
            config_file = 'fibration_invented.ini'
            logging.warning('No .ini file specified --- load %s instead.',
                            config_file)
        config.read(os.path.join(main_file_dir, config_file))
        logging.debug('Configuration sections: %s', config.sections())
        logging.debug('g2 = %s', config.get('fibration', 'g2'))
        logging.debug('g3 = %s', config.get('fibration', 'g3'))

        fibration = EllipticFibration(
            config.get('fibration', 'g2'),
            config.get('fibration', 'g3')
            # config.get('charge', 'fixed_charges'),
            # config.get('charge', 'dsz_matrix'),
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
            config.get('intersection search', 'bin_size')
        )
        kwn.grow(
            # config.get('ODE', 'primary_k_wall_odeint_range'),
            fibration.primary_k_wall_odeint_range,
            config.get('ODE', 'odeint_range'),
            config.get('ODE', 'trajectory_singularity_threshold'),
            config.get('ODE', 'pf_odeint_mxstep'),
            config.get('K-wall network', 'n_iterations'),
            config.get('KSWCF', 'filtration_degree')
        )

        end_time = time.time()
        logging.info('end time: %s', end_time)
        logging.info('elapsed time: %s', end_time - start_time)

        if write_to_file:
            # save picture
            file_path = os.path.join(plots_dir, 
                                     'single_network_' + date_time + '.png')
            save_k_wall_network_plot(kwn, file_path,
                                     plot_range=config.get('plotting', 'range'))
            # save kwn data
            file_path = os.path.join(current_dir, 'results', 
                            'single_network_' + date_time + '.mose')
            saved = f_save(kwn, file_path, config.get('file IO', 
                            'pickle_protocol'))
            print saved
        if show_graphics:
            plot_k_wall_network(
                kwn, config.get('plotting', 'range'),
                plot_bins=show_bins, plot_data_points=show_data_points,
                plot_segments=show_segments
            )

    elif generate_multiple_networks is True:
        if(len(theta_range) == 0):
            theta_range = config.get('MS wall', 'theta_range')
            logging.debug('theta_range = %s from config.', theta_range)
        k_wall_networks = construct_k_wall_networks(
            fibration,
            config.get('intersection search', 'bin_size'),
            # config.get('ODE', 'primary_k_wall_odeint_range'),
            fibration.primary_k_wall_odeint_range,
            config.get('ODE', 'odeint_range'),
            config.get('ODE', 'trajectory_singularity_threshold'),
            config.get('ODE', 'pf_odeint_mxstep'),
            config.get('K-wall network', 'n_iterations'),
            config.get('KSWCF', 'filtration_degree'),
            theta_range,
            config.get('multiprocessing', 'n_processes'),
        )
        ms_walls = build_ms_walls(k_wall_networks)

        end_time = time.time()
        logging.info('end cpu time: %.8f', end_time)
        logging.info('elapsed cpu time: %.8f', end_time - start_time)
             
        if write_to_file:
            # save pictures
            file_path_part = os.path.join(plots_dir, 'phase_scan_' + date_time)
            save_phase_scan(k_wall_networks, ms_walls, file_path_part,
                            plot_range=config.get('plotting', 'range'))
            # save all data
            file_path = os.path.join(current_dir, 'results', 
                            'phase_scan_' + date_time + '.mose')
            saved = f_save([k_wall_networks, ms_walls], file_path,
                           config.get('file IO', 'pickle_protocol'))
            print saved
        if show_graphics:
            plot_ms_walls(ms_walls, plot_range=config.get('plotting', 'range'))


