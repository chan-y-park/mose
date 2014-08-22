from numerics import *

import sympy as sym
from sympy import *
from sympy import mpmath as mp

u = sym.Symbol('u')

g2 = 1 + 4 * u ** 2 / 3
g3 = u / 3 + 8 * u ** 3 / 27

theta = 0

# the options are for numerical integration: initial time, final time, number of steps
options = [0.0, 1.0, 100]

singularties = find_singularities(g2, g3)
print "The Coulomb branch singularities: %s " % singularties

primary_kwalls = grow_primaries(singularties, g2, g3, theta, options)
print "There are %d primary kwalls" % len(primary_kwalls)



