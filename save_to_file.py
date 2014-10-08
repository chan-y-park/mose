from plotting import prepare_k_wall_network_plot, prepare_ms_plot
from matplotlib import pyplot
import pickle
import os


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


def save_k_wall_network_plot(k_wall_network, file_name, plot_range=[],
                             plot_bin=False,
                             plot_intersections=False,
                             display_data_points=False,
                             display_segments=False):
    figure = prepare_k_wall_network_plot(k_wall_network, plot_range, plot_bin,
                                         plot_intersections,
                                         display_data_points,
                                         display_segments)
    pyplot.savefig(file_name)
    return None


def save_ms_plot(ms_walls, plot_range, file_path):
    figure = prepare_ms_plot(ms_walls, plot_range)
    pyplot.savefig(file_path)
    return None


def save_phase_scan(kw_networks, ms_walls, file_path_part, plot_range):
    # First save all k_wall_network snapshots
    for i, k_wall_network in enumerate(kw_networks):
        file_path = file_path_part + '_' + str(i) + '.png'
        save_k_wall_network_plot(k_wall_network, file_path)

    # Then save the plot of ms_walls
    file_path = file_path_part + '_ms_walls.png'
    save_ms_plot(ms_walls, plot_range, file_path)
    return None