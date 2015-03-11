import logging
import os
import sys
import getopt
import time
import pickle
import pdb

import k_wall_network

from matplotlib import pyplot

from elliptic_fibration import EllipticFibration
from k_wall_network import KWallNetwork, construct_k_wall_networks
from k_wall import KWall
from marginal_stability_wall import build_ms_walls


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

    fibration = EllipticFibration(
        config['fibration']['g2'],
        config['fibration']['g3'],
        config['fibration parameters'],
        config['charge']['fixed_charges'],
        config['charge']['dsz_matrix'],
    )

    #KWall.count = 0
    #KWallNetwork.count = 0

    if data['single_network'] is True:
        kwn = KWallNetwork(phase, fibration, config,)
        kwn.grow(config)
        data['k_wall_networks'] = [kwn]

        end_time = time.time()
        logging.info('end time: %s', end_time)
        logging.info('elapsed time: %s', end_time - start_time)

    elif data['multiple_networks'] is True:
        k_wall_networks = construct_k_wall_networks(
            fibration, config,
        )
        data['k_wall_networks'] = k_wall_networks
        ms_walls = build_ms_walls(k_wall_networks)
        data['ms_walls'] = ms_walls

        end_time = time.time()
        logging.info('end cpu time: %.8f', end_time)
        logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return data


def load_data(data_dir):
    """
    Load data from a file
    """
    config = MoseConfig()
    config_file = os.path.join(data_dir, 'config.ini')
    config.read(config_file)
    with open(data_file, 'rb') as fp:
        data = pickle.load(fp)

    return (config, data)


def save_config_file(config, file_name='config.ini', file_dir=None):
    config_file_name = os.path.join(file_dir, file_name)
    logging.info('Save configuration to {}.'.format(config_file_name))
    with open(config_file_name, 'wb') as fp:
        config.parser.write(fp)


def save_data(config, data, data_dir=None, data_file_name='data'):
    data_file_path = os.path.join(data_dir, data_file_name + '.mose')
    with open(data_file_path, 'wb') as fp:
        pickle.dump(data, fp, config['file IO']['pickle_protocol'])
    logging.info('Data saved as' + data_file_path + '.')
