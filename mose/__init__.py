import pdb
import os, platform
import matplotlib
if platform.system() == 'Linux':
    try:
        os.environ['DISPLAY']
        matplotlib.use('TkAgg')
    except KeyError:
        matplotlib.use('Agg') 
else:
    print('Use default backend defined in matplotlibrc: '
          '{}'.format(matplotlib.rcParams['backend']))

from __main__ import run
from api import set_logging, analysis, load, load_config, make_plots, save
set_logging('info')

__all__ = [
    'analysis',
    'load',
    'load_config',
    'make_plots',
    'save',
    'run',
]
