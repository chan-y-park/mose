from plotting import prepare_k_wall_network_plot, prepare_ms_plot
from matplotlib import pyplot
import pickle
import os
# Import necessary modules to read data from saved files
import k_wall_network


def prepare_folder(label):
    """
    prepares the folder where to save data
    """
    # the current directory path
    current_dir = os.getcwd()
    # the general 'results' directory path
    results_dir = os.path.join(current_dir, 'results')
    # If the directory does not exist, create it
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        ### should fix the following
        ###
        # # first, make sure I am not running from INSIDE
        # # the 'mose' folder: I want results to be OUTSIDE
        # if not os.path.exists(os.path.join(current_dir, '__main__.py')):
        #     os.makedirs(results_dir)
        # else:
        #     # this is when we run with ipython from within 
        #     # the mose folder itself
        #     par_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
        #     os.makedirs(os.path.join(par_dir, results_dir))
    # The subdirectory of the plots
    plots_dir = os.path.join(results_dir, label + '_plots')
    os.makedirs(plots_dir)

    return plots_dir


def f_save(obj, file_path, pickle_protocol):
    """
    Save data to a file
    """
    w_stream = open(file_path, 'wb')
    pickle.dump(obj, w_stream, pickle_protocol)
    w_stream.close()
    return "\nContent and picture(s) saved to files " + file_path + "/.png"


def f_recover(file_path):
    """
    Load data from a file
    """
    r_stream = open(file_path, 'rb')
    thawed = pickle.load(r_stream)
    r_stream.close()
    return thawed


def save_k_wall_network_plot(k_wall_network, file_name,
                             plot_bin=False,
                             plot_intersections=False,
                             display_data_points=False,
                             display_segments=False,
                             plot_range=[[-5,5],[-5,5]]
                             ):
    figure = prepare_k_wall_network_plot(k_wall_network, plot_bins=plot_bin,
                                         plot_intersections=plot_intersections,
                                         plot_data_points=display_data_points,
                                         plot_segments=display_segments,
                                         plot_range=plot_range)
    pyplot.savefig(file_name)
    return None


def save_ms_plot(ms_walls, file_path, plot_range=[[-5,5],[-5,5]]):
    figure = prepare_ms_plot(ms_walls, plot_range)
    pyplot.savefig(file_path)
    return None


def save_phase_scan(kw_networks, 
                    ms_walls, 
                    file_path_part, 
                    plot_range=[[-5,5],[-5,5]]
                    ):
    # First save all k_wall_network snapshots
    for i, k_wall_network in enumerate(kw_networks):
        file_path = file_path_part + '_' + str(i) + '.png'
        save_k_wall_network_plot(k_wall_network, file_path,
                                 plot_range=plot_range)

    # Then save the plot of ms_walls
    file_path = file_path_part + '_ms_walls.png'
    save_ms_plot(ms_walls, file_path, plot_range=plot_range)
    return None
