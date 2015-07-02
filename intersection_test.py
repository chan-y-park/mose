import mose
import numpy
import pdb
from mose.plotting import plot_k_walls
from mose.intersection_point import (
    find_new_intersection_points_using_interpolation,
    find_new_intersection_points_using_cgal,
)

def cgal(k_wall_1, k_wall_2):
    return find_new_intersection_points_using_cgal(
        [], [k_wall_1, k_wall_2], [], [[0, -1], [1, 0]], #'Ubuntu'
    )

def find_all_intersections(kws):
    for i, kwi in enumerate(kws):
        for j, kwj in enumerate(kws[i+1:]):
            #plot_k_walls([kwi, kwj])
            print('##################################################')
            print('Finding intersections between K-walls #{} and #{}.'
                  .format(i, i+1+j))
            print('Finding intersections using interpolation:')
            ips = find_new_intersection_points_using_interpolation(
                [], [kwi, kwj], [], [[0, -1], [1, 0]]
            )
            for ip in ips:
                print('intersection at: {}'.format(ip.locus))
            
            ips = cgal(kwi, kwj)
            for ip in ips:
                print ip.locus


#c = mose.load_config('config/fibration_su_2_Nf_1.ini')
#d = mose.analysis(c, phase=1.0)
#d = mose.analysis(c)
#p = mose.make_plots(c, d)
#mose.save(c, d)

c, d = mose.load('results/cgal_test')
kwn = d['k_wall_networks'][0]
kws = kwn.k_walls

find_all_intersections(kws)

#ips1 = cgal(kws[0], kws[2])
#
#for ip in ips1:
#    print ip.locus
#
#ips2 = cgal(kws[0], kws[6])
#
#for ip in ips2:
#    print ip.locus
#
#p = mose.make_plots(c, d)
#kw1 = kws[0]
#kw2 = kws[3]
#ips = mose.intersection_point.find_new_intersection_points(
#    [], [kw1, kw2], [], kwn.fibration.dsz_matrix
#)
#print("intersection at: {}".format(ips[0].locus))
#kws = d['ms_walls'][31].points[0].parents
#ips = mose.intersection_point.find_new_intersection_points(
#    [], kws, [], [[0, -1],[1, 0]]
#)
#msws = d['ms_walls']
#ips = []
#for msw in msws:
#    for point in msw.points:
#        ips.append(point.locus)

#nips = numpy.array(ips)

#with open('k_wall_coordinates.txt', 'w') as f:
#    f.write('K-wall 1:\n')
#    for x, y in kw1.coordinates:
#        f.write('{}, {}\n'.format(x, y))
#    f.write('========================================\n')
#    f.write('K-wall 2:\n')
#    for x, y in kw2.coordinates:
#        f.write('{}, {}\n'.format(x, y))
#

