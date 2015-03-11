import logging
import cmath
import string
import random
import pdb
import sympy
import numpy

#from sympy import Poly
from branch import BranchPoint

NEGLIGIBLE_BOUND = 0.1**12

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class EllipticFibration:
    def __init__(self, g2, g3, params, branch_point_charges, dsz_matrix):

        self.g2 = g2
        self.g3 = g3
        self.params = params
        self.dsz_matrix = dsz_matrix

        self.sym_g2 = sympy.sympify(self.g2)
        self.num_g2 = self.sym_g2.subs(self.params)
        self.sym_g3 = sympy.sympify(self.g3)
        self.num_g3 = self.sym_g3.subs(self.params)

        branch_point_loci = map(
            complex, find_singularities(self.num_g2, self.num_g3, self.params)
        )
        
        ### Introduce string identifiers to label branch-points.
        ### These will be used when building genealogies of intersection 
        ### points, to compare them and build MS walls accordingly.

        bp_identifiers = [id_generator() 
                          for i in range(len(branch_point_loci))]

        # NEED TO DETERMINE THE MONODROMY FROM ANALYSIS OF ELLIPTIC FIBRATION!
        dummy_monodromy = numpy.identity(len(branch_point_charges[0]))

        self.branch_points = [
            BranchPoint(
                branch_point_loci[i],
                branch_point_charges[i], 
                dummy_monodromy,
                bp_identifiers[i]
            )
            for i in range(len(branch_point_loci))
        ]
#        self.branch_cuts = [
#            BranchCut(bp) 
#            for bp in self.branch_points
#        ]


def find_singularities(g2, g3, params):
    """
    find the singularities on the Coulomb branch
    """
    u = sympy.Symbol('u')

    g2_coeffs = map(complex, sympy.Poly(g2, u).all_coeffs())
    g3_coeffs = map(complex, sympy.Poly(g3, u).all_coeffs())
    
    # Converting from
    #   y^2 = 4 x^3 - g_2 x - g_3
    # to 
    #   y^2 = x^3 + f x + g
    
    f = numpy.poly1d(g2_coeffs, variable='u') * (-1 / 4.0)
    g = numpy.poly1d(g3_coeffs, variable='u') * (-1 / 4.0)
    Delta = 4.0 * f ** 3 + 27.0 * g ** 2

    ### Minor Bug: the polynomial is printed with variable 'x' although I 
    ### declared it to be 'u'
    logging.info('discriminant:\n%s', Delta)

    #Accounting for cancellations of higher order terms in discriminant
    for i, coeff in enumerate(Delta.c):
        if numpy.absolute(coeff) > NEGLIGIBLE_BOUND:
            Delta = numpy.poly1d(Delta.c[i:])
            break 

    disc_points = sorted(Delta.r, cmp=lambda x, y: cmp(x.real, y.real))
    logging.info('singularities:\n%s', disc_points)
    return disc_points

