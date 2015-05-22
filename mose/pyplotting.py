import os
import numpy
import logging
import pdb
import matplotlib
from matplotlib import pyplot

import mpldatacursor

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
