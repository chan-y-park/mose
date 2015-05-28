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
from misc import path_derivative, data_plot, complexify


NEGLIGIBLE_BOUND = 0.1**12

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class EllipticFibration:
    def __init__(self, w_f, w_g, params, branch_point_charges=None):

        # self.f = w_f
        # self.g = w_g
        self.params = params

        self.sym_f = sympy.sympify(w_f)
        self.num_f = self.sym_f.subs(self.params)
        self.sym_g = sympy.sympify(w_g)
        self.num_g = self.sym_g.subs(self.params)

        u = sympy.Symbol('u')
        self.f_coeffs = map(complex, sympy.Poly(self.num_f, u).all_coeffs())
        self.g_coeffs = map(complex, sympy.Poly(self.num_g, u).all_coeffs())

        # We will work with the convention that the DSZ matrix is fixed to be
        # the following. Must keep this attribute of the class, as it will be 
        # used when computing the KSWCF for new Kwalls.
        self.dsz_matrix = [[0, -1], [1, 0]]

        ### The following variable is true as long as all the
        ### branch point monodromies are computed using the 
        ### SAME rotation of the x-plane.
        ### It's important that the rotation be the same, because
        ### a rotation may interchange the positions of, say e_1,e_2
        ### resulting in a change of basis in which we express the 
        ### monodromies. As long as we have the same basis for all
        ### monodromies, we are good.
        x_rotation_consistency_check = True

        rotation_n = 10
        zeta = numpy.exp(2 * numpy.pi * 1j / rotation_n)

        for i in range(rotation_n):
            phase = zeta ** (i)
            
            r_num_f = self.sym_f.subs(self.params) / (phase ** 2)
            r_num_g = self.sym_g.subs(self.params) / (phase ** 3)

            r_f_coeffs = map(complex, sympy.Poly(r_num_f, u).all_coeffs())
            r_g_coeffs = map(complex, sympy.Poly(r_num_g, u).all_coeffs())

            self.w_model = wss.WeierstrassModelWithPaths(r_f_coeffs, \
                                                                r_g_coeffs)

            print "I have rotated the x-plane %s times so far." % i
            print "\nf = %s\n\ng = %s\n" % \
                                                (numpy.poly1d(r_f_coeffs), \
                                                numpy.poly1d(r_g_coeffs) \
                                                )

            branch_point_loci = list(self.w_model.get_D().r)

            branch_point_monodromies = [
                wss.monodromy_at_point_wmodel(i, self.w_model) 
                for i in range(len(branch_point_loci))
            ]
            
            x_rotation_consistency_check = \
                                    self.w_model.x_rotation_consistency_check

            ### If, after the i-th rotation, we do have consistency,
            ### we can eventually exit the cycle.
            if x_rotation_consistency_check == True:
                self.x_rotation = phase
                break
            else:
                if i == rotation_n - 1:
                    raise ValueError('Rotated a full 2pi angle, and couldnt '+\
                                    'compute the braiding!')
                else:
                    print "\n\
                    *************************************\n\
                    Having trouble tracking root braiding\n\
                            Will rotate the x-plane      \n\
                    *************************************\n"
        
        ### Now we update the f,g of the fibration class:
        self.sym_f = sympy.sympify(w_f) / (phase ** 2)
        self.num_f = self.sym_f.subs(self.params)
        self.sym_g = sympy.sympify(w_g) / (phase ** 3)
        self.num_g = self.sym_g.subs(self.params)
        self.f_coeffs = map(complex, sympy.Poly(self.num_f, u).all_coeffs())
        self.g_coeffs = map(complex, sympy.Poly(self.num_g, u).all_coeffs())

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

        
        self.branch_points = [
            BranchPoint(
                locus=branch_point_loci[i],
                positive_charge=branch_point_charges[i], 
                monodromy_matrix=branch_point_monodromies[i],
                identifier=bp_identifiers[i],
                fibration = self,
            )
            for i in range(len(branch_point_loci))
        ]

        ### Now determine the "positive periods" of the holomorphic form
        ### on the pinching (gauge) cycle for each branch point.
        ### This will be used to fix the sign ambiguity on the charges of
        ### primary K walls.    
        for bp in self.branch_points:
            ### It is necessary to grow hair before computing the 
            ### reference period, because the latter will be determined 
            ### by comparison with the period at the tip of the hair.
            ### The hair growth will make use of the pf_matrix of the
            ### elliptic fibration, hence of f, g; therefore this 
            ### should be called only after no more rotations of the 
            ### x-plane are necessary.
            print "\nGrowing hair for branch point %s\n" % bp.count
            bp.grow_hair()
            print "\nDetermine positive period for branch point %s\n" % bp.count
            bp.determine_positive_period(reference_period(bp))

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


        def M10(z): 
            return (\
                18432.0 * (f(z) ** 2) * (f_p(z) ** 2) * g_p(z) \
                - 12.0 * f(z) * (1792.0 * g(z) * (f_p(z) ** 3) \
                - 2560.0 * (g_p(z) ** 3)) \
                + 8192.0 * (f(z) ** 3) \
                * (g_p(z) * f_p_p(z) - f_p(z) * g_p_p(z)) \
                + 432.0 * g(z) \
                * (128.0 * g(z) * g_p(z) * f_p_p(z) \
                + (-4.0 * f_p(z)) \
                * ( 16.0 * (g_p(z) ** 2) + 32.0 * g(z) * g_p_p(z))) \
                ) \
                / (16 * (-64.0 * (f(z) ** 3) - 432.0 * (g(z) ** 2)) \
                * (-48.0 * g(z) * f_p(z) + 32.0 * f(z) * g_p(z))) 
        def M11(z):
            return \
            (192.0 * (f(z) ** 2) * f_p(z) + 864.0 * g(z) * g_p(z)) \
            / (-64.0 * (f(z) ** 3) - 432.0 * (g(z) ** 2)) \
            + (16.0 * f_p(z) * g_p(z) + 48.0 * g(z) * f_p_p(z) \
            - 32.0 * f(z) * g_p_p(z)) \
            / (48.0 * g(z) * f_p(z) - 32.0 * f(z) * g_p(z))
        
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

    transpose_m = monodromy.transpose()

    eigen_syst = LA.eig(transpose_m)
    eigen_vals = eigen_syst[0]
    eigen_vects = eigen_syst[1].transpose()

    eigen_val_0 = int(numpy.rint(eigen_vals[0]))
    eigen_val_1 = int(numpy.rint(eigen_vals[1]))
    eigen_vec_0 = numpy.array(eigen_vects[0])[0]
    eigen_vec_1 = numpy.array(eigen_vects[1])[0]

    min_0 = min(abs(eigen_vec_0))
    max_0 = max(abs(eigen_vec_0))
    min_1 = min(abs(eigen_vec_1))
    max_1 = max(abs(eigen_vec_1))
    
    ### A number may be numerically close to zero
    if min_0 < 10**-8:
        min_0 = 0.0
    if max_0 < 10**-8:
        max_0 = 0.0
    if min_1 < 10**-8:
        min_1 = 0.0
    if max_1 < 10**-8:
        max_1 = 0.0

    ### A vector may be of the type (1,0) or (0,1)
    if min_0 == 0.0:
        min_0 = max_0
    if min_1 == 0.0:
        min_1 = max_1

    ### One or both vectors may be (0,0)
    if max_0 == 0.0 and max_1 == 0.0:
        ### both eigenvectors are (0,0)
        print "\nCannot determine eigencharges, they seem to be (0,0).\n"
    elif max_0 == 0.0 and max_1 != 0.0:
        ### the first eigenvector is (0,0), but the second one is not
        norm_eigen_vec_1 = [int(round(x)) for x in list(eigen_vec_1 / min_1)]
        return norm_eigen_vec_1
    elif max_0 != 0.0 and max_1 == 0.0:
        ### the second eigenvector is (0,0), but the first one is not
        norm_eigen_vec_0 = [int(round(x)) for x in list(eigen_vec_0 / min_0)]
        return norm_eigen_vec_0
    else:
        norm_eigen_vec_0 = [int(round(x)) for x in list(eigen_vec_0 / min_0)]
        norm_eigen_vec_1 = [int(round(x)) for x in list(eigen_vec_1 / min_1)]
        n_norm_eigen_vec_1 = [-1 * int(round(x)) \
                                        for x in list(eigen_vec_1 / min_1)]
        if norm_eigen_vec_0 == norm_eigen_vec_1 \
                                    or norm_eigen_vec_0 == n_norm_eigen_vec_1:
            return norm_eigen_vec_0
        else:
            raise ValueError('Cannot compute monodromy eigenvector for:\
                \n'+str(monodromy)+'\ntwo independent eigenvectors!\n')

    


# def positive_period(n, charge, w_model, fibration): 
#     """ 
#     Detemine the period of the holomorphic one form dx / y along the
#     cycle that pinches at a given discriminant locus (the n-th).
#     Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
#     (0,1) determined at the basepoint used to compute all the monodromies.
#     """
    
#     print "\n******************************************\
#            \nEvolving periods for discriminant locus %s\
#            \n******************************************\n" % n

#     ### Setting up the path of integration. Since PF blows up
#     ### at the discriminant locus (although the period does not)
#     ### we will not get exactly to the locus.
#     u0, u1, u2 = w_model.paths_for_periods[n]
#     # path = [u0 + (u1 - u0) * x for x in numpy.linspace(0.0,1.0,1000)] + \
#     #        [u1 + (u2 - u1) * x for x in numpy.linspace(0.0,1.0,1000)]
#     # nint_range = len(path)

    
#     # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
#     period_data = w_model.compute_initial_periods()
#     eta_0 = period_data[0]
#     eta_prime_0 = period_data[1]
#     beta_0 = period_data[2]
#     beta_prime_0 = period_data[3]

#     # print "\nThe integration path\
#     #        \n--------------------\
#     #        \n%s" % path

#     print "\nu_0, eta_0, beta_0:\n%s\n" % [u0, eta_0, beta_0]

#     def deriv(t, y):
#         u, eta, d_eta = y 

#         matrix = fibration.pf_matrix(u)
#         det_pf = abs(det(matrix))
#         # print "u = %s" % u
#         # print "eta = %s" % eta
#         # print "d_eta = %s" % d_eta
#         # print "Jacobian determinant : %s" % det_pf
#         ###
#         ### MOVE THIS PARAMETER ELSEWHERE !!!
#         ###
#         trajectory_singularity_threshold = 10**3
#         ###
#         if det_pf > trajectory_singularity_threshold:
#             singularity_check = True
#             print "\n**************\nhit the singularity threshold!\n**************\n"

#         ### A confusing point to bear in mind: here we are solving the 
#         ### ode with respect to time t, but d_eta is understood to be 
#         ### (d eta / d u), with its own  appropriate b.c. and so on!
#         u_1 = path_derivative(u, u0, u2)
#         eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
#         d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
#         return  array([u_1, eta_1, d_eta_1])
#         # u_1 = (path[int(math.floor(t+1))] - path[int(math.floor(t))]) / 1.0
#         # eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
#         # d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
#         # return  array([u_1, eta_1, d_eta_1])

#     # singularity_check = False
#     ode = scipy.integrate.ode(deriv)
#     ode.set_integrator("zvode")
#     ###
#     ### DEFINE THIS PARAMETER ELSEWHERE !!!
#     ###
#     dt = 0.01
#     ### How many steps before reaching the discriminant locus it should stop
#     ### integrating via PF
#     # cutoff = 0
#     # t1 = len(path) - 2 - cutoff

    
#     recorded_periods = []
#     recorded_loci = []
#     recorded_d_eta = []

#     ### the pinching-cycle period at the starting point, 
#     ### and its derivative
#     eta_gamma_0 = charge[0] * eta_0 + charge[1] * beta_0
#     eta_prime_gamma_0 = charge[0] * eta_prime_0 + charge[1] * beta_prime_0

#     ### Now we PF-evolve the period
#     y_0 = [u0, eta_gamma_0, eta_prime_gamma_0]
#     ode.set_initial_value(y_0)    
#     # while ode.successful() and ode.t < t1 and singularity_check == False:
#     while ode.successful():
#         u, eta, d_eta = ode.y
#         recorded_periods.append(eta)
#         recorded_loci.append(u)
#         recorded_d_eta.append(d_eta)
#         ###
#         ### DEFINE THESE PARAMETERS ELSEWHERE !!!
#         ###
#         if abs(u - u2) < 0.01 or abs(d_eta) > 10:
#             print "\nInterrupting PF transport of period!\n"
#             break
#         else:
#             # print "time: %s" % ode.t
#             ode.integrate(ode.t + dt)

#     u_f, eta_gamma_f, d_eta_gamma_f = ode.y 
#     print "u_f = %s" % u_f
#     print "eta_f = %s\n" % eta_gamma_f

#     data_plot(recorded_loci,"u along PF path")
#     data_plot(recorded_periods,"eta along PF path")
#     data_plot(recorded_d_eta,"d_eta along PF path")

#     return eta_gamma_f


# # def positive_period(n, charge, w_model, fibration): 
    
# #     return 1.0


#### SHOULD MOVE THE FOLLOWING ROUTINE TO THE BRANCH.PY MODULE
def reference_period(branch_point): 
    """ 
    Detemine the period of the holomorphic one form dx / y along the
    cycle that pinches at a given discriminant locus (the n-th).
    Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
    (0,1) determined at the basepoint used to compute all the monodromies.
    """
    
    print "\n******************************************\
           \nEvolving periods for discriminant locus %s\
           \n******************************************\n" % branch_point.count

    charge = branch_point.charge
    fibration = branch_point.fibration
    w_model = fibration.w_model
    
    ### Setting up the path of integration. Since PF blows up
    ### at the discriminant locus (although the period does not)
    ### we will not get exactly to the locus.
    u0 = w_model.base_point
    u2 = complexify(branch_point.hair.coordinates[-1])
    u1 = 1j * u0.imag + u2.real

        
    # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
    period_data = w_model.compute_initial_periods()
    eta_0 = period_data[0]
    eta_prime_0 = period_data[1]
    beta_0 = period_data[2]
    beta_prime_0 = period_data[3]

    # print "\nThe integration path\
    #        \n--------------------\
    #        \n%s" % path

    # print "\nu_0, eta_0, beta_0:\n%s\n" % [u0, eta_0, beta_0]

    def deriv(t, y):
        u, eta, d_eta = y 

        matrix = fibration.pf_matrix(u)
        det_pf = abs(det(matrix))
        # print "u = %s" % u
        # print "eta = %s" % eta
        # print "d_eta = %s" % d_eta
        # print "Jacobian determinant : %s" % det_pf
        ###
        ### MOVE THIS PARAMETER ELSEWHERE !!!
        ###
        trajectory_singularity_threshold = 10**3
        ###
        if det_pf > trajectory_singularity_threshold:
            singularity_check = True
            print "\n**************\nhit the singularity threshold!\n**************\n"

        ### A confusing point to bear in mind: here we are solving the 
        ### ode with respect to time t, but d_eta is understood to be 
        ### (d eta / d u), with its own  appropriate b.c. and so on!
        u_1 = path_derivative(u, u0, u2)
        eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        return  array([u_1, eta_1, d_eta_1])
        # u_1 = (path[int(math.floor(t+1))] - path[int(math.floor(t))]) / 1.0
        # eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        # d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        # return  array([u_1, eta_1, d_eta_1])

    # singularity_check = False
    ode = scipy.integrate.ode(deriv)
    ode.set_integrator("zvode")
    ###
    ### DEFINE THIS PARAMETER ELSEWHERE !!!
    ###
    dt = 0.01
    ### How many steps before reaching the discriminant locus it should stop
    ### integrating via PF
    # cutoff = 0
    # t1 = len(path) - 2 - cutoff

    
    recorded_periods = []
    recorded_loci = []
    recorded_d_eta = []

    ### the pinching-cycle period at the starting point, 
    ### and its derivative
    eta_gamma_0 = charge[0] * eta_0 + charge[1] * beta_0
    eta_prime_gamma_0 = charge[0] * eta_prime_0 + charge[1] * beta_prime_0

    ### Now we PF-evolve the period
    y_0 = [u0, eta_gamma_0, eta_prime_gamma_0]
    ode.set_initial_value(y_0)    
    # while ode.successful() and ode.t < t1 and singularity_check == False:
    while ode.successful():
        u, eta, d_eta = ode.y
        recorded_periods.append(eta)
        recorded_loci.append(u)
        recorded_d_eta.append(d_eta)
        ###
        ### DEFINE THESE PARAMETERS ELSEWHERE !!!
        ###
        if abs(u - u2) < 0.001 or abs(d_eta) > 10:
            # print "\nInterrupting PF transport of period!\n"
            break
        else:
            # print "time: %s" % ode.t
            ode.integrate(ode.t + dt)

    u_f, eta_gamma_f, d_eta_gamma_f = ode.y 
    # print "u_f = %s" % u_f
    # print "eta_f = %s\n" % eta_gamma_f

    data_plot(recorded_loci,"u along reference path")
    data_plot(recorded_periods,"eta along reference path")
    d_eta_by_d_u = [(recorded_periods[i+1] - recorded_periods[i]) \
                    / (recorded_loci[i+1] - recorded_loci[i]) \
                    for i in range(len(recorded_loci)-1)]
    data_plot(recorded_d_eta,"d_eta/d_u along reference path")

    return eta_gamma_f


# def positive_period(n, charge, w_model, fibration): 
    
#     return 1.0



