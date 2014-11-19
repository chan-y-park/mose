import logging
import os
import getopt
import time

from gui import run_gui
from analysis import analysis
from mose_config_parser import MoseConfigParser

shortopts = 'c:l:s:fwgv'
longopts = ['logging_level=',
            'show-bins',
            'show-data-points',
            'show-segments',
]

def run_with_optlist(optlist):

    config_file = ''
    config = MoseConfigParser()
    generate_single_network = False
    generate_multiple_networks = False
    write_to_file = False
    show_graphics = False
    show_bins = False
    show_data_points = False
    show_segments = False
    gui_mode = False
    phase = 0.0
        
    # Default logging
    log = 'warning'
    #logging_level = logging.WARNING
    #logging_format = '%(message)s'

    #theta_range = []

    if len(optlist) == 0:
        print("""usage: python -m mose [OPTION]

    -v:
        run the graphical user interface.

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

    else:

        for opt, arg in optlist:
            if (opt == '-c' and len(arg) > 0):
                config_file = arg
                main_file_dir, main_file_name = os.path.split(__file__)
                config.read(os.path.join(main_file_dir, config_file))
                #theta_range = config.get('MS wall', 'theta_range')
            elif opt == '-v':
                gui_mode = True
            elif (opt == '-l' or opt == '--logging_level'):
                log = arg
            elif opt == '-s':
                # Generate a single K-wall network at a phase
                phase = float(arg)
                analysis_type = 'single'
            elif opt == '-f':
                # Generate K-wall networks at various phases
                analysis_type = 'full'
            elif opt == '-w':
                # save data to external file
                write_to_file = True
            elif opt == '-g':
                # save data to external file
                show_graphics = True
            elif opt == '--show-bins':
                show_bins = True
            elif opt == '--show-data-points':
                show_data_points = True
            elif opt == '--show-segments':
                show_segments = True
        # End of option setting.

        if gui_mode:
            run_gui()
        else:
            analysis(
                        show_graphics,
                        write_to_file,
                        analysis_type,
                        log,
                        config_file,
                        phase,
                        [], #theta_range,
                        show_bins,
                        show_data_points,
                        show_segments,
                    )    

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
