"""
Plots singularities, K-walls, walls of marginal stability, etc., on the
complex 1-dimensional moduli space.
"""

from matplotlib import pyplot

def prepare_k_wall_network_plot(k_wall_network, plot_range=[], plot_bin=False, 
                        plot_intersections=False, 
                        display_data_points=False,
                        display_segments=False): 

    # Give an identifier to the figure we are goint to produce
    pyplot.figure("kwall_snapshot")
    pyplot.figure("kwall_snapshot").clear()

    # Range on the plane to search for intersections
    if(plot_range == []):
        plot_range = k_wall_network.hit_table.get_search_range()
    [[x_min, x_max], [y_min, y_max]] = plot_range

    # Plot setting.
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    # Draw a lattice of bins for visualization.
    if(plot_bin == True):
        bin_size = k_wall_network.hit_table.get_bin_size()
        xv = x_min
        while xv < x_max:
            xv += bin_size
            pyplot.axvline(x = xv, linewidth=0.5, color='0.75')

        yh = y_min  
        while yh < y_max:
            yh += bin_size
            pyplot.axhline(y = yh, linewidth=0.5, color='0.75')
    # End of drawing the bin lattice.

    # Plot branch points
    for bp in k_wall_network.fibration.branch_points:
        bpx = bp.locus.real 
        bpy = bp.locus.imag 
        pyplot.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                    color='k')
    # End of plotting branch points

    # Plot intersection points
    if(plot_intersections == True):
        for ip in k_wall_network.intersections:
            ipx = ip.locus.real 
            ipy = ip.locus.imag 
            pyplot.plot(ipx, ipy, '+', markeredgewidth=2, markersize=8, 
                        color='k')
    # End of plotting intersection points

    # If we have segments of curves, draw them in different colors.
    if(display_segments == True):
        for bin_key in k_wall_network.hit_table:
            for curve_index in hit_table[bin_key]:
                for t_i, t_f in hit_table[bin_key][curve_index]:
                    seg_xcoords, seg_ycoords = \
                        [list(c) for c in \
                            zip(*kwalls[curve_index].coordinates[t_i: t_f + 1])
                        ] 
                    pyplot.plot(seg_xcoords, seg_ycoords, '-')
                    if(display_data_points == True):
                        pyplot.plot(seg_xcoords, seg_ycoords, 'o', color='b')
    else:
        for k_wall in k_wall_network.k_walls:
            xcoords, ycoords = [list(c) for c in zip(*k_wall.coordinates)] 
            pyplot.plot(xcoords, ycoords, '-', color='b')
        if(display_data_points == True):
            pyplot.plot(xcoords, ycoords, 'o', color='b')
    
    return pyplot.figure("kwall_snapshot")



def plot_k_wall_network(k_wall_network, plot_range=[], plot_bin=False, 
                        plot_intersections=False, 
                        display_data_points=False,
                        display_segments=False):
    figure = prepare_k_wall_network_plot(k_wall_network, plot_range, plot_bin, 
                        plot_intersections, display_data_points, 
                        display_segments)
    pyplot.show(figure)
    return None



def prepare_ms_plot(ms_walls, plot_range): 
    """
    Plots MS walls.
    """
    # Give an identifier to the figure we are goint to produce
    pyplot.figure("ms_walls_plot")

    # Range on the plane to search for intersections
    [[x_min, x_max], [y_min, y_max]] = plot_range

    # Plot setting.
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    count = 0
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    # Plot intersection points
    for wall in ms_walls:
        count += 1
        for ip in wall.points:
            ipx = ip.locus.real
            ipy = ip.locus.imag
            pyplot.plot(ipx, ipy, colors[count % len(colors)]+'o', 
                markersize=4)
    # End of plotting intersection points

    return pyplot.figure("ms_walls_plot")
    


def plot_ms_walls(ms_walls, plot_range):
    figure = prepare_ms_plot(ms_walls, plot_range)
    figure.show()
    return None
