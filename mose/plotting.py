"""
Plots singularities, K-walls, walls of marginal stability, etc., on the
complex 1-dimensional moduli space.
"""
import os
import numpy
import Tkinter as tk
import pdb

import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib import pyplot
from matplotlib.widgets import Button
import mpldatacursor

from math import pi

from network_plot import NetworkPlotBase

class KWallNetworkPlotBase(NetworkPlotBase):
    def draw(
        self,
        k_wall_network, 
        plot_range=[[-5, 5], [-5, 5]], 
        plot_joints=False,
        plot_data_points=False,
    ):
        [[x_min, x_max], [y_min, y_max]] = plot_range
        labels = {'branch_points': [], 'joints': [], 'walls': []}

        branch_points = []
        for i, bp in enumerate(k_wall_network.fibration.branch_points):
            branch_points.append([bp.locus.real, bp.locus.imag])
            labels['branch_points'].append("branch point #{}\nM = {}\
                \ncharge = {}".format(i, bp.monodromy_matrix, bp.charge))

        joints = []
        for i, ip in enumerate(k_wall_network.intersections):
            joints.append([ip.locus.real, ip.locus.imag])
            labels['joints'].append("intersection point #{}".format(i))

        walls = []
        for i, k_wall in enumerate(k_wall_network.k_walls):
            xs = k_wall.get_xs()
            ys = k_wall.get_ys()

            k_wall_label = "K-wall #" + str(i) \
                        + "\nDegeneracy: " + str(k_wall.degeneracy) \
                        + "\nIdentifier: " + str(k_wall.identifier)
            num_segments = len(k_wall.splittings) + 1

            x_segments = numpy.split(xs, k_wall.splittings)
            y_segments = numpy.split(ys, k_wall.splittings)
            
            segments = []
            seg_labels = []
            for j in range(num_segments):
                segments.append((x_segments[j], y_segments[j]))
                seg_labels.append(k_wall_label + 
                    "\nLocal charge: {}".format(k_wall.local_charge[j])
                )
            walls.append(segments)
            labels['walls'].append(seg_labels) 
            # else:
            #     labels['walls'].append(k_wall_label)
            # labels['walls'].append("K-wall #{}".format(i))


        super(KWallNetworkPlotBase, self).draw(
            phase=k_wall_network.phase,
            branch_points=branch_points,
            joints=joints,
            walls=walls,
            labels=labels,
            plot_range=[x_min, x_max, y_min, y_max],
            plot_joints=plot_joints,
            plot_data_points=plot_data_points,
        )


class NetworkPlot(KWallNetworkPlotBase):
    """
    This class implements UIs using matplotlib widgets
    so that it can be backend-independent. 
    
    The content of this class is independent of the parent class.
    It only depends on the grandparent class, which should be
    'NetworkPlotBase'. Therefore this class can inherit any class
    whose parent is 'NetworkPlotBase'; just change the name of the
    parent in the definition of this class.
    """
    def __init__(
        self, 
        title=None,
    ):
        super(NetworkPlot, self).__init__(
            matplotlib_figure=pyplot.figure(title),
        )

        self.axes_button_prev = None
        self.axes_button_next = None
        self.index_text = None

    def save(self, plot_dir, file_prefix=''):
        # TODO: change the current figure to plot_id.
        digits = len(str(len(self.plots)-1))
        for i, axes in enumerate(self.plots):
            self.change_current_plot(i)
            plot_file_path = os.path.join(
                plot_dir, file_prefix + str(i).zfill(digits) + '.png'
            )
            self.figure.savefig(plot_file_path)


    def change_current_plot(self, new_plot_idx):
        super(NetworkPlot, self).change_current_plot(new_plot_idx)
        if self.index_text is not None:
            self.index_text.set_text(
                "{}/{}".format(self.current_plot_idx, len(self.plots)-1)
            )



    def show_prev_plot(self, event):
        super(NetworkPlot, self).show_prev_plot(event)


    def show_next_plot(self, event):
        super(NetworkPlot, self).show_next_plot(event)


    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            button_width = .05
            index_width = .03*len(str(len(self.plots)-1))
            button_height = .05
            button_bottom = .025
            center = .5
            margin = .005

            axes_prev_rect = [center - index_width/2 - margin - button_width,
                              button_bottom, button_width, button_height]
            axes_prev = self.figure.add_axes(axes_prev_rect)
            self.button_prev = Button(axes_prev, '<')
            self.button_prev.on_clicked(self.show_prev_plot)

            self.index_text = self.figure.text(
                center - index_width/2, (button_bottom+button_height)/2, 
                "{}/{}".format(self.current_plot_idx, len(self.plots)-1)
            )

            axes_next_rect = [center + index_width/2 + margin,
                              button_bottom, button_width, button_height]
            axes_next = self.figure.add_axes(axes_next_rect)
            self.button_next = Button(axes_next, '>')
            self.button_next.on_clicked(self.show_next_plot)

        self.figure.show()


class NetworkPlotTk(KWallNetworkPlotBase):
    """
    This class implements UIs using Tkinter. 
    
    The content of this class is independent of the parent class.
    It only depends on the grandparent class, which should be
    'NetworkPlotBase'. Therefore this class can inherit any class
    whose parent is 'NetworkPlotBase'; just change the name of the
    parent in the definition of this class.
    """
    def __init__(self, 
        master=None,
        title=None,
    ):
        super(NetworkPlotTk, self).__init__(
            matplotlib_figure=matplotlib.figure.Figure(),
            plot_joints=plot_joints,
            plot_data_points=plot_data_points,
        )

        if master is None:
            master = tk.Tk()
            master.withdraw()
        self.master = master

        # Create a Toplevel widget, which is a child of GUILoom 
        # and contains plots,
        self.toplevel = tk.Toplevel(master)
        self.toplevel.wm_title(title)

        self.plot_idx_scale = None

        self.plot_idx_entry = None
        self.plot_idx_entry_var = tk.StringVar() 
        self.plot_idx_entry_var.trace('w', self.plot_idx_entry_change)

        self.canvas = FigureCanvas(
            self.figure,
            master=self.toplevel,
            resize_callback=self.canvas_resize_callback
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    

    def scale_action(self, scale_value):
        new_plot_idx = int(scale_value)
        self.update_current_plot(new_plot_idx)
        self.plot_idx_entry_var.set(new_plot_idx)


    def plot_idx_entry_change(self, *args):
        try:
            new_plot_idx = int(self.plot_idx_entry_var.get())

            if new_plot_idx == self.current_plot_idx:
                return None
            elif new_plot_idx < 0:
                new_plot_idx = 0
            elif new_plot_idx > len(self.plots) - 1:
                new_plot_idx = len(self.plots) - 1

            self.plot_idx_scale.set(new_plot_idx)
            self.update_current_plot(new_plot_idx)

        except ValueError:
            pass

        return None

    def update_current_plot(self, new_plot_idx):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()

        self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx
        self.set_data_cursor()
        self.canvas.draw_idle()
        self.canvas.get_tk_widget().focus_set()

        return None


    def canvas_resize_callback(self, event):
        self.set_data_cursor()


    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            tk.Label(
                self.toplevel,
                text='Plot #',
            ).pack(side=tk.LEFT)

            self.plot_idx_entry_var.set(plot_idx)
            self.plot_idx_entry = tk.Entry(
                self.toplevel,
                textvariable=self.plot_idx_entry_var,
                width=len(str(len(self.plots)-1)),
            )
            self.plot_idx_entry.pack(side=tk.LEFT)

            tk.Label(
                self.toplevel,
                text='/{}'.format(len(self.plots)-1),
            ).pack(side=tk.LEFT)

            self.plot_idx_scale = tk.Scale(
                self.toplevel,
                command=self.scale_action,
                #length=100*len(self.plots),
                orient=tk.HORIZONTAL,
                showvalue=0,
                to=len(self.plots)-1,
                variable=self.current_plot_idx,
            ) 
            self.plot_idx_scale.pack(
                expand=True,
                fill=tk.X,
                side=tk.LEFT,
            )


class MSWallPlot:
    def __init__(self):
        self.figure = pyplot.figure("Walls of marginal stability") 

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
