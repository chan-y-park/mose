"""
Main module. Run 'python -m mose' at the directory that contains the
directory 'mose'.
"""
import logging
import os
import getopt
import time
import sys
import pdb

from config import MoseConfig
from gui import open_gui
from api import analysis, save_data, load_data, save_config_file
from plotting import KWallNetworkPlot, MSWallPlot
from misc import formatted_date_time


shortopts = 'c:g:hl:p:'
longopts = [
    # options with short versions
    'config-file=',
    'gui-mode=',
    'help',
    'logging_level=',
    'phase=',
    # options without short versions
    'load-data=',
    'save-data=',
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
        'load-data': None,
        'save-data': True,
        'show-plot': False,
    }
        
    for opt, arg in optlist:
        if (opt == '-c' or opt == '--config-file'):
            opts['config-file'] = arg
        elif (opt == '-g' or opt == '--gui-mode'):
            opts['gui-mode'] = eval(arg) 
        elif (opt == '-l' or opt == '--logging_level'):
            opts['logging-level'] = arg
        elif (opt == '-p' or opt == '--phase'):
            # Generate a single K-wall network at a phase
            opts['phase'] = float(arg)
        elif opt == '--load-data':
            opts['load-data'] = arg
        elif opt == '--save-data':
            opts['save-data'] = eval(arg) 
        elif opt == '--show-plot':
            opts['show-plot'] = True
    # End of option setting.

    # Set logging.
    if opts['logging-level'] == 'debug':
        logging_level = logging.DEBUG
        logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
    elif opts['logging-level'] == 'info':
        logging_level = logging.INFO
        logging_format = '%(process)d: %(message)s'
    else:
        logging_level = logging.WARNING
        logging_format = '%(message)s'

    logging.basicConfig(level=logging_level, format=logging_format, 
                        stream=sys.stdout)
    # Load config & data according to the command line args. 
    if opts['load-data'] is not None:
        data_dir = opts['load-data']
        config, data = load_data(data_dir)
    else:
        # Read in the specified config file.
        config = MoseConfig()
        config_file = opts['config-file']
        # No config file chosen; read the default config file.
        if config_file is None:
            logging.warning('No configuration file specified --- '
                            'load "default.ini" instead.')
            config_file = os.path.join(os.curdir, 'default.ini')
        config.read(config_file)
        data = None 

    if opts['gui-mode'] is True:
        return open_gui(config, data)

    if data is None:
        data = analysis(config, phase=opts['phase'],)

    k_wall_network_plot = KWallNetworkPlot()
    ms_wall_plot = MSWallPlot()

    if (opts['show-plot'] or opts['save-data']) is True:
        # Draw the plots of K-wall networks.
        for k_wall_network in data['k_wall_networks']:
            k_wall_network_plot.draw(
                k_wall_network,
                plot_range=config['plotting']['range'], 
            )

        if data['multiple_networks'] is True:
            # Draw MS walls and save the plot.
            ms_wall_plot.draw(
                data['ms_wall'],
                plot_range=config['plotting']['range'], 
            )

    if opts['show-plot'] is True:
        k_wall_network_plot.show()
        if data['multiple_networks'] is True:
            ms_wall_network_plot.show()


    if opts['save-data'] is True:
        # Prepares the folder where to save data.
        current_dir = os.getcwd()
        data_dir = os.path.join(current_dir, 'data', formatted_date_time())
        # If the directory does not exist, create it
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        save_config_file(config, file_dir=data_dir)
        save_data(config, data, data_dir=data_dir)
        k_wall_network_plot.save(data_dir)
        if data['multiple_networks'] is True:
            ms_wall_plot.save(data_dir)

    return (data, k_wall_network_plot, ms_wall_plot)

# Set options from sys.argv when running on the command line,
# then start running the main code.
def run_with_sys_argv(argv):    
    try:
        optlist, args = getopt.getopt(argv, shortopts, longopts,)
        run_with_optlist(optlist)

    except getopt.GetoptError:
        print 'Unknown options.'

# Set options from string 'optstr' when running on the interpreter, 
# then start running the main code.
def run_with_optstr(optstr):
    run_with_sys_argv(optstr.split())


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

--load-data=DATA_FOLDER:
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



if __name__ == '__main__':
    # Set options when running from a command prompt
    run_with_sys_argv(sys.argv[1:])
