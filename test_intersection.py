# This is a test module for intersection.py

import numpy as np
import matplotlib.pyplot as plt
from intersection import *
from math import pi

# Curve 1
# Parameters of the first curve.
x1_min = -2     # minimum value of x range
x1_max = 2      # maximum value of x range
x1_num = 2000   # number of data points
# End of parameters
# Change this functional form to consider different curves.
f1 = lambda x: 1.2*np.cos(2.3*x)

x1 = np.linspace(x1_min, x1_max, x1_num)
y1 = f1(x1)
curve1 = np.transpose([x1, y1]) 
# End of Curve 1

# Cuve 2
x2 = np.linspace(-2, 2, num=3000)
y2 = 0.8 * np.sin(9.3*x2)
curve2 = np.transpose([x2, y2]) 
# End of Curve 2

# Curve 3
# This is a parametric curve.
# Parameters
t3_min = 0
t3_max = 2*pi
t3_num = 2000
# End of parameters
t3 = np.linspace(t3_min, t3_max, t3_num)
x3 = 0.7 * np.cos(t3)
y3 = 0.7 * np.sin(2*t3)
curve3 = np.transpose([x3, y3]) 
# End of Curve 3

# Include your curve in this list to find its intersection with other curves.
lc =[curve1, curve2, curve3] 

# Range on the plane to search for intersections
[[x_min, x_max], [y_min, y_max]] = plane_range = [[-2, 2], [-2, 2]]
# Size of a bin of the coarse-grained plane. Each bin is a square of that size.
bin_size = 0.3

# Draw a lattice of bins for visualization.
xv = x_min
while xv < x_max:
    xv += bin_size
    plt.axvline(x = xv, linewidth=0.5, color='0.75')

yh = y_min  
while yh < y_max:
    yh += bin_size
    plt.axhline(y = yh, linewidth=0.5, color='0.75')
# End of drawing the bin lattice.

# Construct a HitTable to record segmentation of curves.
#ht = HitTable(plane_range, bin_size, bin_fill_offset=5)
ht = HitTable(plane_range, bin_size)
ht.fill(lc)

# Now we have segments of curves. Draw them in different colors.
for bin_key in ht:
    for curve_index in ht[bin_key]:
        for t_i, t_f in ht[bin_key][curve_index]:
            seg_x = [x for x, y in lc[curve_index][t_i: t_f + 1]]   
            seg_y = [y for x, y in lc[curve_index][t_i: t_f + 1]]   
            plt.plot(seg_x, seg_y, '-')

# From the segments find intersection points.
lip = find_intersection(lc, ht, plane_range, bin_size)

# Print and draw intersection points for illustration.
print 'List of intersection points:'
for x, y in lip:
    print '\t[{:+8f}, {:+8f}]'.format(x, y)

for intersection_x, intersection_y in lip:
    plt.plot(intersection_x, intersection_y, '+', markeredgewidth=2, 
                markersize=10, color='k')

# Plot setting.
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.axes().set_aspect('equal')

plt.show()
