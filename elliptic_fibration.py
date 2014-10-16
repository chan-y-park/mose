import logging
import sympy as sym
import numpy as np
import cmath
from sympy import Poly
from branch import BranchPoint, BranchCut

NEGLIGIBLE_BOUND = 0.1**12

class EllipticFibration:
    def __init__(self, g2, g3, branch_point_charges, dsz_matrix,
                    branch_cut_phase, branch_cut_cutoff):

        self.g2 = g2
        self.g3 = g3
        self.dsz_matrix = dsz_matrix
        branch_point_loci = map(complex, find_singularities(g2, g3))
        self.branch_points = [
            BranchPoint(branch_point_loci[i], branch_point_charges[i])  
            for i in range(len(branch_point_loci))
        ]
        self.branch_cuts = [
            BranchCut(bp, branch_cut_phase, branch_cut_cutoff) 
            for bp in self.branch_points
        ]

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

    g2_coeffs = Poly(g2, u).all_coeffs()
    g3_coeffs = Poly(g3, u).all_coeffs()
    
    # Converting from
    # y^2 = 4 x^3 - g_2 x - g_3
    # to 
    # y^2 = x^3 + f x + g
    
    f = np.poly1d(g2_coeffs, variable='u') * (- 1 / 4.0)
    g = np.poly1d(g3_coeffs, variable='u') * (- 1 / 4.0)
    Delta = 4 * f ** 3 + 27 * g ** 2

    ### Minor Bug: the polynomial is printed with variable 'x' although I 
    ### declared it to be 'u'
    logging.info('discriminant: %s', Delta)

    #Accounting for cancellations of higher order terms in discriminant
    for i, coeff in enumerate(Delta.c):
        if np.absolute(coeff) > NEGLIGIBLE_BOUND:
            Delta = np.poly1d(Delta.c[i:])
            break
    
    disc_points = sorted(Delta.r, \
                        cmp=lambda x,y: cmp(x.real, y.real) 
                        )
    logging.info('singularities: %s', disc_points)
    return disc_points

