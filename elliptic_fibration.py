import logging
import sympy as sym
import numpy as np
import cmath
import string
import random
from sympy import Poly
from branch import BranchPoint

NEGLIGIBLE_BOUND = 0.1**12

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class EllipticFibration:
    def __init__(self, g2, g3, branch_point_charges, dsz_matrix):

        self.g2 = g2
        self.g3 = g3
        self.dsz_matrix = dsz_matrix
        branch_point_loci = map(complex, find_singularities(g2, g3))
        
        ### Introduce string identifiers to label branch-points.
        ### These will be used when building genealogies of intersection 
        ### points, to compare them and build MS walls accordingly.

        bp_identifiers = [id_generator() for i in \
                                            range(len(branch_point_loci))]

        # NEED TO DETERMINE THE MONODROMY FROM ANALYSIS OF ELLIPTIC FIBRATION!
        dummy_monodromy = np.identity(len(branch_point_charges[0]))

        self.branch_points = [
            BranchPoint(
                        branch_point_loci[i],
                        branch_point_charges[i], 
                        dummy_monodromy,
                        bp_identifiers[i]
                        )
            for i in range(len(branch_point_loci))
        ]

        # The initial evolution of primary kwalls is handled with an
        # automatic tuning.
        # The length of the single step is calibrated to be
        # 1/2000 th of the minimum distance between any two discriminant
        # loci.
        # The maximal number of steps is set to 400, although it may
        # be automatically truncated whenever e_1, e_2, e_3 become too
        # hard to distinguish.
        step = minimum_distance(self.branch_points) / 2000.0
        max_n_steps = 400
        self.primary_k_wall_odeint_range = [
                                            0.0,\
                                            step * max_n_steps ,\
                                            max_n_steps
                                            ]

        # self.branch_cuts = [
        #     BranchCut(bp) 
        #     for bp in self.branch_points
        # ]

#### OLD METHOD
# def find_singularities(g2, g3):
#     """
#     find the singularities on the Coulomb branch
#     """
#     u = sym.Symbol('u')
#     discriminant = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
#     logging.info('discriminant: %s', discriminant)
#     disc_points = sym.solve(discriminant, u)
#     # disc_points = Poly(discriminant).all_roots()
#     logging.info('singularities: %s', disc_points)
#     return disc_points

#### NEW METHOD
def find_singularities(g2, g3):
    """
    find the singularities on the Coulomb branch
    """
    u = sym.Symbol('u')

    g2_coeffs = map(complex, Poly(g2, u).all_coeffs())
    g3_coeffs = map(complex, Poly(g3, u).all_coeffs())
    
    # Converting from
    # y^2 = 4 x^3 - g_2 x - g_3
    # to 
    # y^2 = x^3 + f x + g
    
    f = np.poly1d(g2_coeffs, variable='u') * (- 1 / 4.0)
    g = np.poly1d(g3_coeffs, variable='u') * (- 1 / 4.0)
    Delta = 4.0 * f ** 3 + 27.0 * g ** 2

    ### Minor Bug: the polynomial is printed with variable 'x' although I 
    ### declared it to be 'u'
    logging.info('discriminant:\n%s', Delta)

    #Accounting for cancellations of higher order terms in discriminant
    for i, coeff in enumerate(Delta.c):
        if np.absolute(coeff) > NEGLIGIBLE_BOUND:
            Delta = np.poly1d(Delta.c[i:])
            break 

    disc_points = sorted(Delta.r, \
                        cmp=lambda x,y: cmp(x.real, y.real) 
                        )
    logging.info('singularities:\n%s', disc_points)
    return disc_points

def minimum_distance(branch_points):
    loci = [bpt.locus for bpt in branch_points]
    min_dist = abs(loci[0]-loci[1])
    for i, z_i in enumerate(loci):
        for j, z_j in list(enumerate(loci))[i+1:]:
            dist = abs(z_i - z_j)
            if dist < min_dist:
                min_dist = dist
    return min_dist

