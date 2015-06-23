import mose
import numpy
from mose.plotting import plot_k_walls

#c = mose.load_config('config/fibration_su_2_Nf_1.ini')
#d = mose.analysis(c, phase=1.0)
#d = mose.analysis(c)
#p = mose.make_plots(c, d)
#mose.save(c, d)

c, d = mose.load('results/test_set')
kwn = d['k_wall_networks'][0]
kws = kwn.k_walls
#p = mose.make_plots(c, d)
kw1 = kws[0]
kw2 = kws[3]
ips = mose.intersection_point.find_new_intersection_points(
    [], [kw1, kw2], [], kwn.fibration.dsz_matrix
)
print("intersection at: {}".format(ips[0].locus))
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
from mose.intersection_point import find_new_intersection_points_using_cgal \
as cgal

intersections = cgal([], [kw1, kw2], [], [[0, -1], [1, 0]], 'Ubuntu')
for ip in intersections:
    print ip.locus
