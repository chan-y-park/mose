import pdb
import os, platform
import time
import matplotlib
if matplotlib.rcParams['backend'] == 'nbAgg':
    print('Use IPython notebook backend for matplotlib.')
elif platform.system() == 'Linux':
    try:
        os.environ['DISPLAY']
        matplotlib.use('TkAgg')
        print('Use TkAgg backend for matplotlib.')
    except KeyError:
        matplotlib.use('Agg') 
        print('Use Agg backend for matplotlib.')
else:
    print('Use default backend defined in matplotlibrc: '
          '{}'.format(matplotlib.rcParams['backend']))

from api import set_logging, load_config, make_plots
from api import generate_k_wall_network as generate
from api import load_k_wall_network as load
from api import save_k_wall_network as save

MOSE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..',)

LOGGING_FILE_PATH = os.path.join(
    MOSE_DIR,
    ('logs/mose_{}-{:02}-{:02} {:02}:{:02}:{:02}.log'
     .format(*time.localtime(time.time())[:6])),
)

set_logging('info', logging_file_name=LOGGING_FILE_PATH)
#set_logging('debug', logging_file_name=LOGGING_FILE_PATH)

__all__ = [
    'load_config',
    'generate'
    'make_plots',
    'load',
    'save',
]
