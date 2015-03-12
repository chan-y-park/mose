"""
Plots singularities, K-walls, walls of marginal stability, etc., on the
complex 1-dimensional moduli space.
"""
import os
import numpy
import matplotlib
import logging
import Tkinter as tk
import pdb
import matplotlib
# use() directive must be called before importing matplotlib.pyplot
matplotlib.use('TkAgg')     

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib import pyplot


class KWallNetworkPlot:
    def __init__(self, root=None):
        # Create a Toplevel widget, which is a child of GUILoom 
        # and contains plots,
        self.root = root
        if root is None:
            self.toplevel = tk.Tk()
        else:
            self.toplevel = tk.Toplevel(root)
        self.toplevel.wm_title('K-wall Network Plot')
        self.plots = []
        self.current_plot_idx = None 

        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(
            self.figure,
            master=self.toplevel,
            #resize_callback=(
            #    lambda event: self.set_data_cursor()
            #)
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack()

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack()

    def draw(
        self,
        k_wall_network, 
        plot_range=[[-5, 5], [-5, 5]], 
        plot_bins=False, 
        plot_intersections=False, 
        plot_data_points=False,
        plot_segments=False,
    ):
        hit_table = k_wall_network.hit_table
        k_walls = k_wall_network.k_walls
        phase = k_wall_network.phase

        [[x_min, x_max], [y_min, y_max]] = plot_range

        rect = [0.125, 0.15, 0.8, 0.75]
        axes = self.figure.add_axes(
            rect,
            label='kwn_snapshot@{}.'.format(phase),    
            xlim=(x_min, x_max),
            ylim=(y_min, y_max),
            aspect='equal',
        )

        # Draw a lattice of bins for visualization.
        if(plot_bins is True):
            bin_size = hit_table.get_bin_size()
            xv = x_min
            while xv < x_max:
                xv += bin_size
                axes.axvline(x = xv, linewidth=0.5, color='0.75')

            yh = y_min  
            while yh < y_max:
                yh += bin_size
                axes.axhline(y = yh, linewidth=0.5, color='0.75')
        # End of drawing the bin lattice.

        # Plot branch points
        for bp in k_wall_network.fibration.branch_points:
            bpx = bp.locus.real 
            bpy = bp.locus.imag 
            axes.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                        color='k')
        # End of plotting branch points

        # Plot intersection points
        if(plot_intersections == True):
            for ip in k_wall_network.intersections:
                ipx = ip.locus.real 
                ipy = ip.locus.imag 
                axes.plot(ipx, ipy, '+', markeredgewidth=2, markersize=8, 
                            color='k')
        # End of plotting intersection points

        # If we have segments of curves, draw them in different colors.
        if(plot_segments == True):
            for bin_key in hit_table:
                for curve_index in hit_table[bin_key]:
                    for t_i, t_f in hit_table[bin_key][curve_index]:
                        seg_xcoords, seg_ycoords = [
                            list(c) for c in
                            zip(*k_walls[curve_index].coordinates[t_i: t_f+1])
                        ] 
                        axes.plot(seg_xcoords, seg_ycoords, '-')
                        if(plot_data_points == True):
                            axes.plot(seg_xcoords, seg_ycoords, 'o',
                                        color='b')
        else:
            for k_wall in k_walls:
                xcoords, ycoords = [list(c) for c in zip(*k_wall.coordinates)] 
                axes.plot(xcoords, ycoords, '-', color='b')
            if(plot_data_points == True):
                axes.plot(xcoords, ycoords, 'o', color='b')
        
        self.plots.append(axes)
        if self.current_plot_idx is not None:
            self.plots[self.current_plot_idx].set_visible(False)
            self.current_plot_idx += 1
        else:
            self.current_plot_idx = 0
        axes.set_visible(True)

        return None


    def save(self, plot_dir, file_prefix='k_wall_network_'):
        # TODO: change the current figure to plot_id.
        digits = len(str(len(self.plots)-1))
        for i, axes in enumerate(self.plots):
            plot_file_path = os.path.join(
                plot_dir, file_prefix + str(i).zfill(digits) + '.png'
            )
            self.figure.savefig(plot_file_path)

    def show(self):
        self.plots[self.current_plot_idx].set_visible(True)
        if self.root is None:
            self.toplevel.mainloop()


class MSWallPlot:
    def __init__(self):
        self.figure = matplotlib.figure.Figure() 

    def draw(
        self,
        ms_walls,
        plot_range=[[-5,5],[-5,5]]
    ):
        """
        Plots MS walls.
        """
        # Range on the plane to search for intersections
        [[x_min, x_max], [y_min, y_max]] = plot_range

        rect = [0.125, 0.15, 0.8, 0.75]
        axes = self.figure.add_axes(
            rect,
            label='ms_walls',    
            xlim=(x_min, x_max),
            ylim=(y_min, y_max),
            aspect='equal',
        )

        count = 0
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

        # Plot intersection points
        for wall in ms_walls:
            count += 1
            for ip in wall.points:
                ipx = ip.locus.real
                ipy = ip.locus.imag
                axes.plot(ipx, ipy, colors[count % len(colors)]+'o', 
                    markersize=4)

    def show(self):
        self.figure.show()


    def save(self, plot_dir, file_prefix='ms_walls'):
        plot_file_path = os.path.join(
            plot_dir, file_prefix + '.png'
        )
        self.figure.savefig(plot_file_path)


def plot_coordinates(coordinates, plot_range=[[-5, 5], [-5, 5]],
                     plot_data_points=False,):
    # Plot setting.
    [[x_min, x_max], [y_min, y_max]] = plot_range
    pyplot.figure()
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    xcoords, ycoords = [list(c) for c in zip(*coordinates)] 
    pyplot.plot(xcoords, ycoords, '-', color='b')

    if(plot_data_points == True):
        pyplot.plot(xcoords, ycoords, 'o', color='k', markersize=4)

    pyplot.show()

def plot_eta(eta):
    pyplot.figure()
    X = numpy.arange(-10, 10, 0.1)
    Y = numpy.arange(-10, 10, 0.1)
    X, Y = numpy.meshgrid(X, Y)
    Z_min = -5.0
    Z_max = 5.0
    n_levels = 20
    levels = [Z_min + (Z_max - Z_min)*i/n_levels for i in range(n_levels)]

    eta_r = numpy.vectorize(lambda z: complex(eta(z)).real)
    eta_i = numpy.vectorize(lambda z: complex(eta(z)).imag)
    Z_r = eta_r(X + 1j * Y)
    Z_i = eta_i(X + 1j * Y)

    pyplot.subplot(1, 2, 1, aspect=1)
    csr = pyplot.contourf(X, Y, Z_r, levels=levels)
    pyplot.colorbar(csr, shrink=.5)

    pyplot.subplot(1, 2, 2, aspect=1)
    csi = pyplot.contourf(X, Y, Z_i, levels=levels)
    pyplot.colorbar(csi, shrink=.5)

    pyplot.show()
