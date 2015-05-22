import numpy
import pdb
import logging
import mpldatacursor

import matplotlib
#matplotlib.use('Agg')     

from math import pi

class AggNetworkPlot(object):
    def __init__(self, 
        matplotlib_figure=None,
        plot_joints=False,
        plot_data_points=False,
    ):
        self.plot_joints = plot_joints
        self.plot_data_points = plot_data_points

        self.plots = []
        self.data_cursor = None
        self.current_plot_idx = None 
        self.axes_button_prev = None
        self.axes_button_next = None
        self.index_text = None

        self.figure = matplotlib_figure
    

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
        rect = [.1, 0.15, .8, .8]
        #rect = [0.125, 0.125, .8, 0.8]

        axes = self.figure.add_axes(
            rect,
            label="Network #{}".format(len(self.plots)),
            xlim=(x_min, x_max),
            ylim=(y_min, y_max),
            aspect='equal',
        )

        axes.set_title('phase = ({:.4f})pi'.format(phase/pi))

        for i, wall in enumerate(walls):
            xs = wall.get_xs()
            ys = wall.get_ys()
            if(self.plot_data_points is True):
                axes.plot(xs, ys, 'o', color='k')
            else:
                axes.plot(xs, ys, '-',
                          #color='b',
                          label=labels['walls'][i],)
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

        axes.set_visible(False)
        self.plots.append(axes)


    def set_data_cursor(self):
        if self.current_plot_idx is None:
            return None

        # Use a DataCursor to interactively display the label
        # for artists of the current axes.
        self.data_cursor = mpldatacursor.datacursor(
            axes=self.plots[self.current_plot_idx],
            formatter='{label}'.format,
            tolerance=4,
            hover=True,
            #display='single',
            #display='multiple',
            #draggable=True,
        )


    def change_current_plot(self, new_plot_idx):
        if new_plot_idx == self.current_plot_idx:
            return None
        elif new_plot_idx < 0:
            new_plot_idx = 0
        elif new_plot_idx > len(self.plots) - 1:
            new_plot_idx = len(self.plots) - 1

        self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx
        if self.index_text is not None:
            self.index_text.set_text(
                "{}/{}".format(self.current_plot_idx, len(self.plots)-1)
            )


    def show_prev_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx-1)
        self.set_data_cursor()
        self.figure.show()


    def show_next_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx+1)
        self.set_data_cursor()
        self.figure.show()


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

