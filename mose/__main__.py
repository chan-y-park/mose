"""
Main module. Run 'python -m mose' at the directory that contains the
directory 'mose'.
"""
import logging
import os
import getopt
import time
import sys
import Tkinter as tk
import pdb

from config import MoseConfig
from gui import open_gui
from api import (set_logging, analysis, save, load, load_config, make_plots,
                run_diagnostics)


shortopts = 'c:g:hl:p:'
longopts = [
    # options with short versions
    'config-file=',
    'gui-mode=',
    'help',
    'logging_level=',
    'phase=',
    # options without short versions
    'load-data',
    'load-data-from=',
    'save-data',
    #'save-data-at=',
    'show-plot',
]


def run_with_optlist(optlist):
    if len(optlist) == 0:
        return print_help()

    opts = {
        'config-file': None,
        'gui-mode': False,
        'logging-level': 'info',
        'phase': None,
        'load-data': False,
        'load-data-from': None,
        'save-data': False,
        #'save-data-at': None,
        'show-plot': False,
    }

    for opt, arg in optlist:
        if (opt == '-h' or opt == '--help'):
            return print_help()
        elif (opt == '-c' or opt == '--config-file'):
            opts['config-file'] = arg
        elif (opt == '-g' or opt == '--gui-mode'):
            opts['gui-mode'] = eval(arg) 
        elif (opt == '-l' or opt == '--logging_level'):
            opts['logging-level'] = arg
        elif (opt == '-p' or opt == '--phase'):
            # Generate a single K-wall network at a phase
            opts['phase'] = float(arg)
        elif opt == '--load-data':
            opts['load-data'] = True 
        elif opt == '--load-data-from':
            opts['load-data-from'] = arg
        elif opt == '--save-data':
            opts['save-data'] = True 
        #elif opt == '--save-data-at':
        #    opts['save-data-at'] = arg 
        elif opt == '--show-plot':
            opts['show-plot'] = True
    # End of option setting.

    set_logging(opts['logging-level'])

    if opts['gui-mode'] is True:
        return open_gui(config, data)

    # Load config & data. 
    if opts['load-data'] is True:
        # Load config & data using a file dialog. 
        config, data = load()    
        if data is None:
            return None
    elif opts['load-data-from'] is not None:
        # Load config & data according to the command line args. 
        data_dir = opts['load-data-from']
        config, data = load(data_dir)
    else:
        # Read in the specified config file.
        config = load_config(opts['config-file'])
        if config is None:
            return None
        data = None 

    if data is None:
        data = analysis(config, phase=opts['phase'],)

    if (
        opts['show-plot'] is True
        # (opts['show-plot'] or opts['save-data']) is True
        #or (opts['save-data-at'] is not None)
    ):
        k_wall_network_plot, ms_wall_plot = make_plots(
            config, data, show_plot=opts['show-plot'],
        )

    if opts['save-data'] is True:
        #save(config, data, k_wall_network_plot, ms_wall_plot)
        save(config, data, open_dialog=False)


def run(optstr='', argv=None):
    # Set options from string 'optstr' when running on the interpreter, 
    # then start running the main code.
    if argv is None:
        argv = optstr.split()

    # Set options from sys.argv (when running on the command line)
    # or from argv from optstr (when running on the interpreter),
    # then start running the main code.
    try:
        optlist, args = getopt.getopt(argv, shortopts, longopts,)
        return run_with_optlist(optlist)

    except getopt.GetoptError as e:
        print 'Unknown option: {}.'.format(e)
        return None


def print_help():
    print("""usage: python -m mose [OPTION]

-c CFG_FILE_NAME, --config-file=CFG_FILE_NAME:
    read CFG_FILE_NAME to set up the configuration.
    The file path is relative to the execution directory.

-g, --gui-mode:
    run the graphical user interface.

-l LEVEL, --logging_level=LEVEL:
    set logging level to LEVEL.

-p THETA, --phase=THETA:
    produce a single K-wall network at the phase of THETA
    and plot the K-wall network.

--load-data:
    open a dialog window to choose a folder containing
    configuration & data files to load.

--load-data-from=DATA_FOLDER:
    load data from a folder containing configuration & data
    files.

--save-data:
    save data on a date-named file for offline analysis,
    the name will include 'phase_scan' or 'single_network'
    according to what data is stored.
    Will also produce saved pictures.
    
--show-plot:
    Show the plots of the generated data or the loaded data.
    For data at a single phase, it will show a plot of 
    the K-wall network; for data at multiple phases, it
    will also show the plot of marginal stability walls.
    """)

    return None



if __name__ == '__main__':
    # Set options when running from a command prompt
    run(sys.argv[1:])
