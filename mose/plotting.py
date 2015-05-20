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
import mpldatacursor

from matplotlib import pyplot
from network_plot import NetworkPlot


class KWallNetworkPlot(NetworkPlot):
    def __init__(self, matplotlib_figure=None, plot_joints=False, 
                 plot_data_points=False):
        super(KWallNetworkPlot, self).__init__(
            matplotlib_figure=matplotlib_figure, plot_joints=plot_joints,
            plot_data_points=plot_data_points,
        )

    def draw(
        self,
        k_wall_network, 
        plot_range=[[-5, 5], [-5, 5]], 
    ):
        [[x_min, x_max], [y_min, y_max]] = plot_range
        branch_points = []
        joints = []
        labels = {'branch_points': [], 'joints': [], 'walls': []}
        for i, bp in enumerate(k_wall_network.fibration.branch_points):
            branch_points.append([bp.locus.real, bp.locus.imag])
            labels['branch_points'].append("branch point #{}".format(i))
        for i, ip in enumerate(k_wall_network.intersections):
            joints.append([ip.locus.real, ip.locus.imag])
            labels['joints'].append("intersection point #{}".format(i))
        for i, wall in enumerate(k_wall_network.k_walls):
            kwall_label = "K-wall #" + str(i) \
                        + "\nInitial charge: " + str(wall.charge(0)) \
                        + "\nDegeneracy: " + str(wall.degeneracy)
            labels['walls'].append(kwall_label)
            # labels['walls'].append("K-wall #{}".format(i))

        super(KWallNetworkPlot, self).draw(
            phase=k_wall_network.phase,
            branch_points=branch_points,
            joints=joints,
            walls=k_wall_network.k_walls,
            labels=labels,
            plot_range=[x_min, x_max, y_min, y_max],
        )


    def save(self, plot_dir, file_prefix='k_wall_network_'):
        # TODO: change the current figure to plot_id.
        digits = len(str(len(self.plots)-1))
        for i, axes in enumerate(self.plots):
            self.change_current_plot(i)
            plot_file_path = os.path.join(
                plot_dir, file_prefix + str(i).zfill(digits) + '.png'
            )
            self.figure.savefig(plot_file_path)


class MSWallPlot:
    def __init__(self, matplotlib_figure=None):
        self.figure = matplotlib_figure

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

        #count = 0
        #colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

        # Plot intersection points
        for i, wall in enumerate(ms_walls):
            label = "MS wall #{}".format(i)
            xs = []
            ys = []
            for l in wall.locus:
                xs.append(l.real)
                ys.append(l.imag)
            axes.plot(xs, ys, '-', markersize=4, label=label)
            axes.plot(xs, ys, 'o', markersize=4, label=label)

        data_cursor = mpldatacursor.datacursor(
            axes=axes,
            formatter='{label}'.format,
            #tolerance=2,
            hover=True,
            #display='single',
            #display='multiple',
        )


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
