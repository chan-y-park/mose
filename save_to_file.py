from config import PICKLE_PROTOCOL
from plotting import prepare_k_wall_network_plot, prepare_ms_plot
from matplotlib import pyplot
import pickle



def f_save(obj, file_name):
    """
    Save data to a file
    """
    w_stream = open(file_name, 'wb')
    pickle.dump(obj, w_stream, PICKLE_PROTOCOL)
    w_stream.close()
    return "\nContent and picture(s) saved to files " + file_name + "/.png"


def f_recover(file_name):
    """
    Load data from a file
    """
    r_stream = open(file_name,'rb')
    thawed = pickle.load(r_stream)
    r_stream.close()
    return thawed


def save_k_wall_network_plot(k_wall_network, file_name, plot_range=[],
                        plot_bin=False, 
                        plot_intersections=False, 
                        display_data_points=False,
                        display_segments=False):
    figure = prepare_k_wall_network_plot(k_wall_network, plot_range, plot_bin, 
                        plot_intersections, display_data_points, 
                        display_segments)
    pyplot.savefig(file_name)
    return None


def save_ms_plot(ms_walls, plot_range, file_name):
    figure = prepare_ms_plot(ms_walls, plot_range)
    pyplot.savefig(file_name)
    return None


def save_phase_scan(kw_networks, ms_walls, file_name_part, plot_range):
    ### First save all k_wall_network snapshots
    for i, k_wall_network in enumerate(kw_networks):
        file_name = file_name_part + '_' + str(i) + '.png'
        save_k_wall_network_plot(k_wall_network, file_name)

    ### Then save the plot of ms_walls
    file_name = file_name_part + '_ms_walls.png'
    save_ms_plot(ms_walls, plot_range, file_name)
    return None
