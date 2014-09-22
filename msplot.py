"""
Plots singularities, K-walls, walls of marginal stability, etc., on the
complex 1-dimensional moduli space.
"""

from matplotlib import pyplot

def kwallplot(kwalls, plot_range, bin_size=0, intersection_points=[],
                hit_table={}, mark_data_plot=False): 

    # Range on the plane to search for intersections
    [[x_min, x_max], [y_min, y_max]] = plot_range

    # Plot setting.
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    # Draw a lattice of bins for visualization.
    if(bin_size > 0):
        xv = x_min
        while xv < x_max:
            xv += bin_size
            pyplot.axvline(x = xv, linewidth=0.5, color='0.75')

        yh = y_min  
        while yh < y_max:
            yh += bin_size
            pyplot.axhline(y = yh, linewidth=0.5, color='0.75')
    # End of drawing the bin lattice.

    # Plot intersection points
    if(len(intersection_points)>0):
        for ip in intersection_points:
            ipx = ip.locus.real 
            ipy = ip.locus.imag 
            pyplot.plot(ipx, ipy, '+', markeredgewidth=2, markersize=10, 
                        color='k')
    # End of plotting intersection points

    # If we have segments of curves, draw them in different colors.
    if(len(hit_table)>0):
        for bin_key in hit_table:
            for curve_index in hit_table[bin_key]:
                for t_i, t_f in hit_table[bin_key][curve_index]:
                    seg_xcoords, seg_ycoords = \
                        [list(c) for c in \
                            zip(*kwalls[curve_index].coordinates[t_i: t_f + 1])
                        ] 
                    pyplot.plot(seg_xcoords, seg_ycoords, '-')
                    if(mark_data_plot == True):
                        pyplot.plot(seg_xcoords, seg_ycoords, 'o', color='b')
    else:
        for trajectory in kwalls:
            xcoords, ycoords = [list(c) for c in zip(*trajectory.coordinates)] 
            pyplot.plot(xcoords, ycoords, '-', color='b')
        if(mark_data_plot == True):
            pyplot.plot(xcoords, ycoords, 'o', color='b')

    pyplot.show()



