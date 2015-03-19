import logging
import math
import cmath
import string
import random
import pdb
import sympy
import numpy

import weierstrass as wss

#from sympy import Poly
from branch import BranchPoint

NEGLIGIBLE_BOUND = 0.1**12

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class EllipticFibration:
    def __init__(self, g2, g3, params):

        self.g2 = g2
        self.g3 = g3
        self.params = params

        self.sym_g2 = sympy.sympify(self.g2)
        self.num_g2 = self.sym_g2.subs(self.params)
        self.sym_g3 = sympy.sympify(self.g3)
        self.num_g3 = self.sym_g3.subs(self.params)

        u = sympy.Symbol('u')
        g2_coeffs = map(complex, sympy.Poly(self.num_g2, u).all_coeffs())
        g3_coeffs = map(complex, sympy.Poly(self.num_g3, u).all_coeffs())
        # Converting from
        #   y^2 = 4 x^3 - g_2 x - g_3
        # to 
        #   y^2 = x^3 + f x + g
        f_coeffs = numpy.poly1d(g2_coeffs, variable='u') * (-1 / 4.0)
        g_coeffs = numpy.poly1d(g3_coeffs, variable='u') * (-1 / 4.0)

        self.w_model = wss.WeierstrassModelWithPaths(f_coeffs, g_coeffs)

        # We will work with the convention that the DSZ matrix is fixed to be
        # the following. Must keep this attribute of the class, as it will be 
        # used when computing the KSWCF for new Kwalls.
        self.dsz_matrix = [[0, 1], [-1, 0]]

        
        ### RE-ENABLE THIS ###
        # branch_point_loci = list(self.w_model.get_D())
        branch_point_loci = map(
            complex, find_singularities(self.num_g2, self.num_g3, self.params)
        )

        ### RESTORE THIS!
        # branch_point_monodromies = \
        #                 [wss.monodromy_at_point_wmodel(i, self.w_model) \
        #                 for i in range(len(branch_point_loci))]
        # THIS IS A TEMPORARY DUMMY ASSIGNMENT
        branch_point_monodromies = [
            [[1, -2],[0, 1]] for i in range(len(branch_point_loci))
        ]
        
        branch_point_charges = [
            monodromy_eigencharge(m) for m in branch_point_monodromies
        ]

        
        ### Introduce string identifiers to label branch-points.
        ### These will be used when building genealogies of intersection 
        ### points, to compare them and build MS walls accordingly.

        bp_identifiers = [id_generator() 
                          for i in range(len(branch_point_loci))]

        print """
              Here is the full list of branch point data:
              -------------------------------------------
              """
        for i in range(len(branch_point_loci)):
            print (
                """
                Branch point # {}
                --------------------
                Locus: {}
                Monodromy: {}
                Charge: {}
                Internal label: {}
                """.format(
                    i,
                    branch_point_loci[i],
                    branch_point_monodromies[i],
                    branch_point_charges[i],
                    bp_identifiers[i],
                )
            )
        ### Now determine the "positive periods" of the holomorphic form
        ### on the pinching (gauge) cycle for each branch point.
        ### This will be used to fix the sign ambiguity on the charges of
        ### primary K walls.            
        branch_point_positive_periods = [
            positive_period(i, branch_point_charges[i], self.w_model, self)
            for i in range(len(branch_point_loci))
        ]


        self.branch_points = [
            BranchPoint(
                locus=branch_point_loci[i],
                positive_charge=branch_point_charges[i], 
                monodromy_matrix=branch_point_monodromies[i],
                positive_period=branch_point_positive_periods[i],
                identifier=bp_identifiers[i],
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
            0.0, step * max_n_steps, max_n_steps
        ]

#        self.branch_cuts = [
#            BranchCut(bp) 
#            for bp in self.branch_points
#        ]


#### NOW SUPERSEDED BY WEIERSTRASS CLASS ITSELF
def find_singularities(num_g2, num_g3, params):
    """
    find the singularities on the Coulomb branch
    """
    u = sympy.Symbol('u')

    g2_coeffs = map(complex, sympy.Poly(num_g2, u).all_coeffs())
    g3_coeffs = map(complex, sympy.Poly(num_g3, u).all_coeffs())
    
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

def minimum_distance(branch_points):
    loci = [bpt.locus for bpt in branch_points]
    min_dist = abs(loci[0]-loci[1])
    for i, z_i in enumerate(loci):
        for j, z_j in list(enumerate(loci))[i+1:]:
            dist = abs(z_i - z_j)
            if dist < min_dist:
                min_dist = dist
    return min_dist

def monodromy_eigencharge(monodromy):
    # m = magnetic charge
    # n = electric charge
    # we work in conventions of Seiberg-Witten 2, with monodromy matrices
    # acting from the LEFT (hence our matrices are TRANSPOSE os those in SW!), 
    # and given by:
    #       1 + 2 n m    &    - m^2       \\
    #       4 n^2        &    1 - 2 n m

    nm = (monodromy[0][0] - monodromy[1][1]) / 4.0
    m_temp = math.sqrt(-1.0 * monodromy[0][1])
    n_temp = 2 * math.sqrt(monodromy[1][0] / 4.0)
    if nm != 0:
        if nm > 0:
            n = n_temp
            m = m_temp
        else:
            n = -n_temp
            m = m_temp
    else:
        m = m_temp
        n = n_temp
    return (int(m), int(n))


# def positive_period(n, charge, w_model, fibration): 
#     """ 
#     Detemine the period of the holomorphic one form dx / y along the
#     cycle that pinches at a given discriminant locus (the n-th).
#     Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
#     (0,1) determined at the basepoint used to compute all the monodromies.
#     """
    
#     ### NEED TO EXTRACT THIS DATA FROM W-MODEL
#     path = ...
#     nint_range = len(path)
#     u_0 = path[0]
#     # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
#     eta_0 = ...
#     eta_prime_0 = ...
#     beta_0 = ...
#     beta_prime_0 = ...

#     # Setting up generic data for PF evolution
#     g2 = fibration.g2
#     g3 = fibration.g3
#     u = sympy.Symbol('u')
#     x = sympy.Symbol('x')
#     Delta = sympy.simplify(g2 ** 3 - 27 * g3 ** 2)
#     delta = sympy.simplify(3 * (g3) * diff(g2, u) - 2 * (g2) * diff(g3, u))

#     pf_matrix = lambdify(u, 
#         [
#             [0, 1], 
#             [
#                 sympy.simplify(
#                     (-(3*g2*(delta**2) / (16*(Delta**2))) + 
#                     ((diff(delta,u) * diff(Delta,u)) / (12*delta*Delta)) - 
#                     ((diff(Delta, u, 2)) / (12*Delta)) + 
#                     (((diff(Delta, u))**2) / (144 * (Delta**2) ) ) ) 
#                 ), 
#                 sympy.simplify( 
#                     ( (diff(delta,u) / delta) - (diff(Delta,u) / Delta) ) 
#                 ) 
#             ]
#         ] 
#     )
    
#     def deriv(y,t):
#         z = y[0] + 1j * y[1]
#         eta = y[2] + 1j * y[3]
#         d_eta = y[4] + 1j * y[5]
#         matrix = pf_matrix(z)
#         det_pf = abs(det(matrix))
#         # print "PF matrix = %s" % matrix
#         # print "PF determinant = %s" % det_pf
#         if det_pf > trajectory_singularity_threshold:
#             self.singular = True
#             self.singular_point = z

#         # A confusing point to bear in mind: here we are solving the 
#         # ode with respect to time t, but d_eta is understood to be 
#         # (d eta / d u), with its own  appropriate b.c. and so on!

#         # Note here I am adapting the derivative of z to be precisely the
#         # one determined by the path on the Coulomb branch.
#         z_1 = (path[int(math.floor(t))] - path[int(math.floor(t+1))]) / 1.0
#         eta_1 = z_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
#         d_eta_1 = z_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
#         return  array(
#                 [
#                     z_1.real, z_1.imag, 
#                     eta_1.real, eta_1.imag, 
#                     d_eta_1.real, d_eta_1.imag
#                 ]
#         )
#         # the following rescaled matrix will not blow up at singularities, 
#         # but it will not reach the singularities either..
#         # return abs(det(matrix)) * array([z_1.real, z_1.imag, 
#         #              eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])


#     # Now first we PF-evolve the period of (1,0)
#     singularity_check= False
#     t0 = 0
#     y0 = map(N, 
#         array(
#             [
#                 (u_0).real,
#                 (u_0).imag,
#                 (eta_0).real,
#                 (eta_0).imag,
#                 (eta_prime_0).real,
#                 (eta_prime_0).imag
#             ]
#         )    
#     )
#     time = linspace(*nint_range)
#     y = odeint(deriv, y0, time)
#     eta_final = y[2] + 1j * y[3]

#     # Second we PF-evolve the period of (0,1)
#     singularity_check= False
#     t0 = 0
#     y0 = map(N, 
#         array(
#             [
#                 (u_0).real,
#                 (u_0).imag,
#                 (beta_0).real,
#                 (beta_0).imag,
#                 (beta_prime_0).real,
#                 (beta_prime_0).imag
#             ]
#         )    
#     )
#     time = linspace(*nint_range)
#     y = odeint(deriv, y0, time)
#     beta_final = y[2] + 1j * y[3]    

#     positive_period_final = charge[0] * eta_final + charge[1] * beta_final
#     return positive_period_final


def positive_period(n, charge, w_model, fibration): 
    
    return 1.0
