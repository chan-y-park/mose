import numpy
import pdb
import logging
import Tkinter as tk
import mpldatacursor

import matplotlib
# use() directive must be called before importing matplotlib.pyplot
matplotlib.use('TkAgg')     

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib import pyplot
from math import pi

class TkNetworkPlot(object):
    def __init__(self, 
        master=None,
        #plot_on_cylinder=False,
        #plot_bins=False, 
        plot_joints=False,
        plot_data_points=False,
        #plot_segments=False,
    ):
        if master is None:
            master = tk.Tk()
            master.withdraw()

        self.master = master
        #self.plot_on_cylinder = plot_on_cylinder
        #self.plot_bins = plot_bins
        self.plot_joints = plot_joints
        self.plot_data_points = plot_data_points
        #self.plot_segments = plot_segments

        # Create a Toplevel widget, which is a child of GUILoom 
        # and contains plots,
        self.toplevel = tk.Toplevel(master)
        #self.toplevel.wm_title('Spectral Network Plot')

        self.plots = []
        self.data_cursor = None
        self.current_plot_idx = None 

        self.plot_idx_scale = None

        self.plot_idx_entry = None
        self.plot_idx_entry_var = tk.StringVar() 
        self.plot_idx_entry_var.trace('w', self.plot_idx_entry_change)

        self.figure = matplotlib.figure.Figure()
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
    

    def draw(self, phase=None, branch_points=None, joints=None, walls=None,
             labels=None, plot_range=[-5, 5, -5, 5]):
        """
        branch_points = [[bpx, bpy], ...]
        joints = [[jpx, jpy], ...]
        walls = [[wall.get_xs(), wall.get_ys(), ...]
        labels = {'branch_points': [bp1_label, ...],
                  'joints': [jp1_label, ...],
                  'walls': [wall1_label, ...]}
        """
        x_min, x_max, y_min, y_max = plot_range
        rect = [0.125, 0.15, .8, 0.75]

        axes = self.figure.add_axes(
            rect,
            label=phase,
            xlim=(x_min, x_max),
            ylim=(y_min, y_max),
            aspect='equal',
        )

        axes.set_title('phase = ({:.4f})pi'.format(phase/pi))

        # Plot branch points
        for i, bp in enumerate(branch_points):
            bpx, bpy = bp
            axes.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                      color='k', label=labels['branch_points'][i],)
        # End of plotting branch points
   
        # Plot joints
        if(self.plot_joints is True):
            for i, jp in enumerate(joints):
                jpx, jpy = jp
                axes.plot(jpx, jpy, '+', markeredgewidth=2,
                          markersize=8, color='k', label=labels['joints'][i],)
        # End of plotting joints

        for i, wall in enumerate(walls):
            xs = wall.get_xs()
            ys = wall.get_ys()
            if(self.plot_data_points is True):
                axes.plot(xs, ys, 'o', color='k')
            else:
                axes.plot(xs, ys, '-',
                          #color='b',
                          label=labels['walls'][i],)
        axes.set_visible(False)
        self.plots.append(axes)

        return None


    def set_data_cursor(self):
        if self.current_plot_idx is None:
            return None

        # Use a DataCursor to interactively display the label
        # for artists of the current axes.
        self.data_cursor = mpldatacursor.datacursor(
            axes=self.plots[self.current_plot_idx],
            formatter='{label}'.format,
            tolerance=2,
            #hover=True,
            #display='single',
            display='multiple',
        )
    
        return None 


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
            self.data_cursor.hide()

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
