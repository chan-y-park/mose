"""
Plots singularities, K-walls, walls of marginal stability, etc., on the
complex 1-dimensional moduli space.
"""
import os
import numpy
import pdb

from agg_network_plot import AggNetworkPlot
from tk_network_plot import TkNetworkPlot

use_tk = True


class KWallNetworkPlot(TkNetworkPlot if use_tk is True else AggNetworkPlot):
    def __init__(
        self, master=None, plot_joints=False, plot_data_points=False,
    ):
        super(KWallNetworkPlot, self).__init__(
            master=master,
            title='K-wall network plot',
            plot_joints=plot_joints,
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
            labels['branch_points'].append("branch point #{}\nM = {}\
                \ncharge = {}".format(i, bp.monodromy_matrix, bp.charge))
        for i, ip in enumerate(k_wall_network.intersections):
            joints.append([ip.locus.real, ip.locus.imag])
            labels['joints'].append("intersection point #{}".format(i))
        for i, k_wall in enumerate(k_wall_network.k_walls):
            k_wall_label = "K-wall #" + str(i) \
                        + "\nDegeneracy: " + str(k_wall.degeneracy) \
                        + "\nIdentifier: " + str(k_wall.identifier)
            num_segments = len(k_wall.splittings)

            if num_segments == 0:
                single_label = k_wall_label + \
                    "\nLocal charge: {}".format(k_wall.initial_charge)
                labels['walls'].append(single_label) 
            elif num_segments > 0:
                seg_labels = []
                for j in range(num_segments+1):
                    seg_labels.append(k_wall_label + 
                        "\nLocal charge: {}".format(k_wall.local_charge[j])
                    )
                labels['walls'].append(seg_labels) 
            # else:
            #     labels['walls'].append(k_wall_label)
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

