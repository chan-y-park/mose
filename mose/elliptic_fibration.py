import logging
import math
import cmath
import string
import random
import pdb
import sympy
import numpy
import scipy

import weierstrass as wss

from cmath import pi, exp, phase, sqrt
from numpy import linalg as LA

from branch import BranchPoint

from scipy.integrate import odeint
from scipy import interpolate

# temporarily added the following:
from sympy import diff, N, simplify
from sympy import mpmath as mp
from sympy.utilities.lambdify import lambdify
from cmath import exp, pi
from numpy import array, linspace
from numpy.linalg import det
from operator import itemgetter


NEGLIGIBLE_BOUND = 0.1**12

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class EllipticFibration:
    def __init__(self, w_f, w_g, params, branch_point_charges=None):

        self.f = w_f
        self.g = w_g
        self.params = params

        self.sym_f = sympy.sympify(self.f)
        self.num_f = self.sym_f.subs(self.params)
        self.sym_g = sympy.sympify(self.g)
        self.num_g = self.sym_g.subs(self.params)

        u = sympy.Symbol('u')
        self.f_coeffs = map(complex, sympy.Poly(self.num_f, u).all_coeffs())
        self.g_coeffs = map(complex, sympy.Poly(self.num_g, u).all_coeffs())
        
        # Converting from
        #   y^2 = 4 x^3 - g_2 x - g_3
        # to 
        #   y^2 = x^3 + f x + g
        # f_coeffs = numpy.poly1d(self.g2_coeffs, variable='u') * (-1 / 4.0)
        # g_coeffs = numpy.poly1d(self.g3_coeffs, variable='u') * (-1 / 4.0)

        self.w_model = wss.WeierstrassModelWithPaths(self.f_coeffs, \
                                                     self.g_coeffs)

        # We will work with the convention that the DSZ matrix is fixed to be
        # the following. Must keep this attribute of the class, as it will be 
        # used when computing the KSWCF for new Kwalls.
        self.dsz_matrix = [[0, 1], [-1, 0]]
       
        branch_point_loci = list(self.w_model.get_D().r)

        
        ### Calculate branch point charges using weierstrss.py
        ### If trouble is encountered, can resort to the patch commented below.
        ###
        branch_point_monodromies = [
            wss.monodromy_at_point_wmodel(i, self.w_model) 
            for i in range(len(branch_point_loci))
        ]
        
        branch_point_charges = [
            monodromy_eigencharge(m) for m in branch_point_monodromies
        ]

        ### The following might be useful if some trouble is encountered
        ### with automatic computation of discriminant charges from weierstrass 
        ### module.
        ###
        # if branch_point_charges is None:
        #     # Calculate branch point charges using weierstrss.py
        #     branch_point_monodromies = [
        #         wss.monodromy_at_point_wmodel(i, self.w_model) 
        #         for i in range(len(branch_point_loci))
        #     ]
            
        #     branch_point_charges = [
        #         monodromy_eigencharge(m) for m in branch_point_monodromies
        #     ]
        # ### ONCE THE WEIERSTRASS MODULE WORKS RELIABLY, NEED TO REMOVE THE 
        # ### 'IF' STRUCTURE, THE FOLLOWING PART WON'T BE NEEDED.
        # else:
        #     # THIS IS A TEMPORARY DUMMY ASSIGNMENT
        #     branch_point_monodromies = [
        #         [[1, -2],[0, 1]] for i in range(len(branch_point_loci))
        #     ]

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
                    numpy.array(branch_point_monodromies[i]).tolist(),
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

#        self.branch_cuts = [
#            BranchCut(bp) 
#            for bp in self.branch_points
#        ]
    
    def pf_matrix(self, z):
        
        f = numpy.poly1d(self.f_coeffs)
        g = numpy.poly1d(self.g_coeffs)
        f_p = f.deriv()
        g_p = g.deriv()
        f_p_p = f_p.deriv()
        g_p_p = g_p.deriv()

        ### Once done with replacing conventions successfully, 
        ### should simply the following expression.

        def M10(z): 
            return (\
                -18 * ((-4.0 * f(z)) ** 2) * ((-4.0 * f_p(z)) ** 2) * (-4.0 * g_p(z)) \
                + 3 * (-4.0 * f(z)) * (7 * (-4.0 * g(z)) * ((-4.0 * f_p(z)) ** 3) \
                + 40 * ((-4.0 * g_p(z)) ** 3)) \
                + ((-4.0 * f(z)) ** 3) * (-8 * (-4.0 * g_p(z)) * (-4.0 * f_p_p(z)) \
                + 8 * (-4.0 * f_p(z)) * (-4.0 * g_p_p(z))) \
                -108 * (-4.0 * g(z)) \
                * (-2 * (-4.0 * g(z)) * (-4.0 * g_p(z)) * (-4.0 * f_p_p(z)) \
                + (-4.0 * f_p(z)) \
                * (((-4.0 * g_p(z)) ** 2) + 2 * (-4.0 * g(z)) * (-4.0 * g_p_p(z)))) \
                ) \
                / (16 * (((-4.0 * f(z)) ** 3) -27 * ((-4.0 * g(z)) ** 2)) \
                * (-3 * (-4.0 * g(z)) * (-4.0 * f_p(z)) + 2 * (-4.0 * f(z)) * (-4.0 * g_p(z)))) 
        def M11(z):
            return \
            (-3 * ((-4.0 * f(z)) ** 2) * (-4.0 * f_p(z)) + 54 * (-4.0 * g(z)) * (-4.0 * g_p(z))) \
            / (((-4.0 * f(z)) ** 3) - (27 * (-4.0 * g(z)) ** 2)) \
            + ((-4.0 * f_p(z)) * (-4.0 * g_p(z)) + 3 * (-4.0 * g(z)) * (-4.0 * f_p_p(z)) \
            - 2 * (-4.0 * f(z)) * (-4.0 * g_p_p(z))) \
            / (3 * (-4.0 * g(z)) * (-4.0 * f_p(z)) - 2 * (-4.0 * f(z)) * (-4.0 * g_p(z)))
        
        return [[0, 1], [M10(z), M11(z)]]


        ### Attempt to use polyval, to make evaluation faster, 
        ### needs some extra work.
        ###
        # print "TYPE: %s" % type((\
        #         -18 * (g2 ** 2) * (g2_p ** 2) * g3_p \
        #         + 3 * g2 * (7 * g3 * (g2_p ** 3) \
        #         + 40 * (g3_p ** 3)) \
        #         + (g2 ** 3) * (-8 * g3_p * g2_p_p \
        #         + 8 * g2_p * g3_p_p) \
        #         -108 * g3 \
        #         * (-2 * g3 * g3_p * g2_p_p \
        #         + g2_p \
        #         * ((g3_p ** 2) + 2 * g3 * g3_p_p)) \
        #         ) \
        #         / (16 * ((g2 ** 3) -27 * (g3 ** 2)) \
        #         * (-3 * g3 * g2_p + 2 * g2 * g3_p)\
        #         ))

        # M10 = numpy.polyval(\
        #         (\
        #         -18 * (g2 ** 2) * (g2_p ** 2) * g3_p \
        #         + 3 * g2 * (7 * g3 * (g2_p ** 3) \
        #         + 40 * (g3_p ** 3)) \
        #         + (g2 ** 3) * (-8 * g3_p * g2_p_p \
        #         + 8 * g2_p * g3_p_p) \
        #         -108 * g3 \
        #         * (-2 * g3 * g3_p * g2_p_p \
        #         + g2_p \
        #         * ((g3_p ** 2) + 2 * g3 * g3_p_p)) \
        #         ) \
        #         / (16 * ((g2 ** 3) -27 * (g3 ** 2)) \
        #         * (-3 * g3 * g2_p + 2 * g2 * g3_p)\
        #         ), z)

        # M11 = numpy.polyval(\
        #     (-3 * (g2 ** 2) * g2_p + 54 * g3 * g3_p) \
        #     / ((g2 ** 3) - (27 * g3 ** 2)) \
        #     + (g2_p * g3_p + 3 * g3 * g2_p_p \
        #     - 2 * g2 * g3_p_p) \
        #     / (3 * g3 * g2_p - 2 * g2 * g3_p), z)

        # print "PF matrix : %s" % [[0, 1], [M10, M11]]
        # return [[0, 1], [M10, M11]]


#### NOW SUPERSEDED BY WEIERSTRASS CLASS ITSELF
#def find_singularities(num_g2, num_g3, params):
#    """
#    find the singularities on the Coulomb branch
#    """
#    u = sympy.Symbol('u')
#
#    g2_coeffs = map(complex, sympy.Poly(num_g2, u).all_coeffs())
#    g3_coeffs = map(complex, sympy.Poly(num_g3, u).all_coeffs())
#    
#    # Converting from
#    #   y^2 = 4 x^3 - g_2 x - g_3
#    # to 
#    #   y^2 = x^3 + f x + g
#    
#    f = numpy.poly1d(g2_coeffs, variable='u') * (-1 / 4.0)
#    g = numpy.poly1d(g3_coeffs, variable='u') * (-1 / 4.0)
#    Delta = 4.0 * f ** 3 + 27.0 * g ** 2
#
#    ### Minor Bug: the polynomial is printed with variable 'x' although I 
#    ### declared it to be 'u'
#    logging.info('discriminant:\n%s', Delta)
#
#    #Accounting for cancellations of higher order terms in discriminant
#    for i, coeff in enumerate(Delta.c):
#        if numpy.absolute(coeff) > NEGLIGIBLE_BOUND:
#            Delta = numpy.poly1d(Delta.c[i:])
#            break 
#
#    disc_points = sorted(Delta.r, cmp=lambda x, y: cmp(x.real, y.real))
#    logging.info('singularities:\n%s', disc_points)
#    return disc_points


def monodromy_eigencharge(monodromy):
    # # m = magnetic charge
    # # n = electric charge
    # # we work in conventions of Seiberg-Witten 2, with monodromy matrices
    # # acting from the LEFT (hence our matrices are TRANSPOSE os those in SW!), 
    # # and given by:
    # #       1 + 2 n m    &    - m^2       \\
    # #       4 n^2        &    1 - 2 n m

    # # print "this is the monodromy matrix: \n%s" % monodromy
    # # print "this is the monodromy matrix type: %s" % type(monodromy)

    # nm = (monodromy[0,0] - monodromy[1,1]) / 4.0
    # m_temp = math.sqrt(-1.0 * monodromy[0,1])
    # n_temp = 2 * math.sqrt(monodromy[1,0] / 4.0)
    # if nm != 0:
    #     if nm > 0:
    #         n = n_temp
    #         m = m_temp
    #     else:
    #         n = -n_temp
    #         m = m_temp
    # else:
    #     m = m_temp
    #     n = n_temp
    # return (int(m), int(n))

    eigen_sys = LA.eig(monodromy)
    eigen_val_0 = int(numpy.rint(eigen_sys[0][0]))
    eigen_val_1 = int(numpy.rint(eigen_sys[0][1]))
    min_0 = min(abs(numpy.array(eigen_sys[1][0])[0]))
    eigen_vec_0 = (numpy.array(eigen_sys[1][0])[0] / min_0).astype(int).tolist()
    min_1 = min(abs(numpy.array(eigen_sys[1][1])[0]))
    eigen_vec_1 = (numpy.array(eigen_sys[1][1])[0] / min_1).astype(int).tolist()
    n_eigen_vec_1 = (-1 * numpy.array(eigen_sys[1][1])[0] / min_1).astype(int).tolist()

    if eigen_vec_0 == eigen_vec_1 or eigen_vec_0 == n_eigen_vec_1:
        return eigen_vec_0
    else:
        raise ValueError('Cannot compute monodromy eigenvector for the matrix\
            \n'+str(monodromy))


def positive_period(n, charge, w_model, fibration): 
    """ 
    Detemine the period of the holomorphic one form dx / y along the
    cycle that pinches at a given discriminant locus (the n-th).
    Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
    (0,1) determined at the basepoint used to compute all the monodromies.
    """
    
    print "\n******************************************\
           \nEvolving periods for discriminant locus %s\
           \n******************************************\n" % n

    ### Setting up the path of integration. Since PF blows up
    ### at the discriminant locus (although the period does not)
    ### we will not get exactly to the locus.
    u_0, u_1, u_2 = w_model.paths_for_periods[n]
    path = [u_0 + (u_1 - u_0) * x for x in numpy.linspace(0.0,1.0,1000)] + \
           [u_1 + (u_2 - u_1) * x for x in numpy.linspace(0.0,1.0,1000)]
    nint_range = len(path)

    
    # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
    period_data = w_model.compute_initial_periods()
    eta_0 = period_data[0]
    eta_prime_0 = period_data[1]
    beta_0 = period_data[2]
    beta_prime_0 = period_data[3]

    # print "\nThe integration path\
    #        \n--------------------\
    #        \n%s" % path

    print "\nu_0, eta_0, beta_0:\n%s\n" % [u_0, eta_0, beta_0]

    def deriv(t, y):
        u, eta, d_eta = y 

        matrix = fibration.pf_matrix(u)
        det_pf = abs(det(matrix))
        ###
        ### MOVE THIS PARAMETER ELSEWHERE !!!
        ###
        trajectory_singularity_threshold = 10**3
        ###
        if det_pf > trajectory_singularity_threshold:
            singularity_check = True

        # A confusing point to bear in mind: here we are solving the 
        # ode with respect to time t, but d_eta is understood to be 
        # (d eta / d u), with its own  appropriate b.c. and so on!
        u_1 = (path[int(math.floor(t+1))] - path[int(math.floor(t))]) / 1.0
        eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        return  array([u_1, eta_1, d_eta_1])


    singularity_check = False
    ode = scipy.integrate.ode(deriv)
    ode.set_integrator("zvode")
    dt = 1
    ### How many steps before reaching the discriminant locus it should stop
    ### integrating via PF
    cutoff = 0
    t1 = len(path) - 2 - cutoff

    ### the pinching-cycle period at the starting point, 
    ### and its derivative
    eta_gamma_0 = charge[0] * eta_0 + charge[1] * beta_0
    eta_prime_gamma_0 = charge[0] * eta_prime_0 + charge[1] * beta_prime_0

    ### Now we PF-evolve the period
    y_0 = [u_0, eta_gamma_0, eta_prime_gamma_0]
    ode.set_initial_value(y_0)    
    while ode.successful() and ode.t < t1 and singularity_check == False:
        # print "time: %s" % ode.t
        ode.integrate(ode.t + dt)

    u_f, eta_gamma_f, d_eta_gamma_f = ode.y 
    print "u_f = %s" % u_f
    print "eta_f = %s\n" % eta_gamma_f

    return eta_gamma_f


# def positive_period(n, charge, w_model, fibration): 
    
#     return 1.0
