#--------------------------------------------------------
### Temporary knobs to try different combinations
### of the available algorithms

ENHANCE_MS_WALLS = True

# If false, will use charge orbits
SORT_BY_GENEALOGY = True

# If false, will use straight approximation
CUT_K_WALLS = False

# Three options: 'sweep', 'neighbor', 'phase'
MS_WALLS_SORTING = 'phase'

# Whether to connect intersection points of an MS wall with a line
PLOT_MS_WALL_LINKS = True

# Decide whether to ignore intersections with pairing >2 or <0
# will make the code run much faster
IGNORE_WILD_INTERSECTIONS = True

#--------------------------------------------------------


import logging
import os
import sys
import getopt
import time
import pickle
import Tkinter as tk
import tkFileDialog
import shutil
import pdb

import matplotlib
import k_wall_network

from config import MoseConfig
from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from k_wall import KWall
from marginal_stability_wall import build_ms_walls
from plotting import NetworkPlot, NetworkPlotTk
from plotting import MSWallPlot
from misc import formatted_date_time
from diagnostics import diagnose_kwall_network

LOGGING_FILE_NAME = 'log.mose.txt'



def set_logging(level):
    if level == 'debug':
        logging_level = logging.DEBUG
        logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
    elif level == 'info':
        logging_level = logging.INFO
        logging_format = '%(process)d: %(message)s'
    else:
        logging_level = logging.WARNING
        logging_format = '%(message)s'

    # logging.basicConfig(level=logging_level, format=logging_format, 
    #                     stream=sys.stdout)
    
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    ### create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    formatter = logging.Formatter(logging_format)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Log to a file 'log.mose.txt'
    fh = logging.FileHandler(LOGGING_FILE_NAME, 'w')
    fh.setLevel(logging_level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)


def analysis(config, phase=None,):
    data = {
        'k_wall_networks': None,
        'ms_walls': None,
        'single_network': None,
        'multiple_networks': None,
    }
    if phase is not None:
        data['single_network'] = True
        data['multiple_networks'] = False
    else:
        data['single_network'] = False 
        data['multiple_networks'] = True

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    trajectory_singularity_threshold = (
            config['ODE']['trajectory_singularity_threshold']
        )
    ode_size_of_step = config['ODE']['size_of_step']
    ode_num_steps = config['ODE']['num_steps']

    fibration = EllipticFibration(
        config['fibration']['f'],
        config['fibration']['g'],
        config['fibration parameters'],
        config['KSWCF']['dsz_matrix']
    )

    #KWall.count = 0
    #KWallNetwork.count = 0

    if data['single_network'] is True:
        kwn = KWallNetwork(phase, fibration, config)
        kwn.grow(config)
        data['k_wall_networks'] = [kwn]

        end_time = time.time()
        logging.info('end time: %s', end_time)
        logging.info('elapsed time: %s', end_time - start_time)

    elif data['multiple_networks'] is True:
        k_wall_networks = construct_k_wall_networks(
            fibration, config
        )
        data['k_wall_networks'] = k_wall_networks
        ms_walls = build_ms_walls(k_wall_networks)
        data['ms_walls'] = ms_walls

        end_time = time.time()
        logging.info('end cpu time: %.8f', end_time)
        logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return data


def analysis_with_parameters(config, phase=None, parameters=None):
    data = {
        'k_wall_networks': None,
        'ms_walls': None,
        'single_network': None,
        'multiple_networks': None,
    }
    if phase is not None:
        data['single_network'] = True
        data['multiple_networks'] = False
    else:
        data['single_network'] = False 
        data['multiple_networks'] = True

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    trajectory_singularity_threshold = (
            config['ODE']['trajectory_singularity_threshold']
        )
    ode_size_of_step = config['ODE']['size_of_step']
    ode_num_steps = config['ODE']['num_steps']

    fibration = EllipticFibration(
        config['fibration']['f'],
        config['fibration']['g'],
        parameters
    )

    #KWall.count = 0
    #KWallNetwork.count = 0

    if data['single_network'] is True:
        kwn = KWallNetwork(phase, fibration, config)
        kwn.grow(config)
        data['k_wall_networks'] = [kwn]

        end_time = time.time()
        logging.info('end time: %s', end_time)
        logging.info('elapsed time: %s', end_time - start_time)

    elif data['multiple_networks'] is True:
        k_wall_networks = construct_k_wall_networks(
            fibration, config
        )
        data['k_wall_networks'] = k_wall_networks
        ms_walls = build_ms_walls(k_wall_networks)
        data['ms_walls'] = ms_walls

        end_time = time.time()
        logging.info('end cpu time: %.8f', end_time)
        logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return data


def run_diagnostics(kwn):
    diagnose_kwall_network(kwn)
    

def load_config(config_file=None):
    config = MoseConfig()
    if config_file is None:
        root = tk.Tk()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': root,
            'title': 'Select a configuration file to load.',
        }
        config_file = tkFileDialog.askopenfilename(**file_opts)
        root.destroy()
        if config_file == '':
            return None
    else:
        config_file = os.path.join(os.curdir, config_file)

    config.read(config_file)
    return config


def load(data_dir=None):
    """
    Load config & data from a given directory
    """
    if data_dir is None:
        root = tk.Tk()
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': root,
            'title': 'Select a directory that contains data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        root.destroy()
        if data_dir == '':
            return (None, None)

    logging.info('Opening data directory "{}"...'.format(data_dir))

    config_file = os.path.join(data_dir, 'config.ini')
    config = load_config(config_file)
    data_file = os.path.join(data_dir, 'data.mose')
    with open(data_file, 'rb') as fp:
        data = pickle.load(fp)

    return (config, data)


def save_config(config, file_name='config.ini', file_dir=None):
    config_file_name = os.path.join(file_dir, file_name)
    logging.info('Save configuration to {}.'.format(config_file_name))
    with open(config_file_name, 'wb') as fp:
        config.parser.write(fp)


def save_data(config, data, data_dir=None, data_file_name='data'):
    data_file_path = os.path.join(data_dir, data_file_name + '.mose')
    with open(data_file_path, 'wb') as fp:
        pickle.dump(data, fp, config['file IO']['pickle_protocol'])
    logging.info("Data saved as {}.".format(data_file_path))


def save(config, data, k_wall_network_plot=None, ms_wall_plot=None,
         data_dir=None, open_dialog=False):
    """
    Save config & data files in a directory.
    """
    if open_dialog is True:
        root = tk.Tk()
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': root,
            'title': 'Select a directory to save files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        root.destroy()
        if data_dir == '':
            return None 

    if data_dir is None:
        # Prepares the folder where to save data.
        current_dir = os.getcwd()
        data_dir = os.path.join(current_dir, 'results', formatted_date_time())
    # If the directory does not exist, create it
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    save_config(config, file_dir=data_dir)
    save_data(config, data, data_dir=data_dir)
    if k_wall_network_plot is not None:
        k_wall_network_plot.save(data_dir, file_prefix='k_wall_network_')
    if ms_wall_plot is not None:
        ms_wall_plot.save(data_dir)
    # Copy log file to data_dir
    shutil.copyfile(LOGGING_FILE_NAME,
                    os.path.join(data_dir, LOGGING_FILE_NAME))


def make_plots(config, data, show_plot=True, master=None,
               k_wall_network_plot_range=None, ms_wall_plot_range=None):
    k_wall_network_plot_title = 'K-wall Network'
    if matplotlib.rcParams['backend'] == 'TkAgg':
        k_wall_network_plot = NetworkPlotTk(
            title=k_wall_network_plot_title,
        )
    else:
        k_wall_network_plot = NetworkPlot(
            title=k_wall_network_plot_title,
        )

    # Draw the plots of K-wall networks.
    for k_wall_network in data['k_wall_networks']:
        k_wall_network_plot.draw(
            k_wall_network,
            #plot_range=config['plotting']['range'], 
            plot_range=k_wall_network_plot_range, 
        )
        if k_wall_network_plot_range is None:
            k_wall_network_plot.autoscale()

    if data['multiple_networks'] is True:
        ms_wall_plot = MSWallPlot()
        # Draw MS walls and save the plot.
        ms_wall_plot.draw(
            data['ms_walls'],
            #plot_range=config['plotting']['range'], 
            plot_range=ms_wall_plot_range, 
        )
    else:
        ms_wall_plot = None

    if show_plot is True:
        k_wall_network_plot.show()
        if ms_wall_plot is not None:
            ms_wall_plot.show()
        try:
            raw_input('Press any key to continue...')
        except NameError:
            pass

#    if master is None:
#        try:
#            raw_input('Press any key to continue...')
#        except NameError:
#            pass
    return (k_wall_network_plot, ms_wall_plot)
