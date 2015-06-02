#                            NUMERICAL CONSTANTS
#------------------------------------------------------------------------
### Defines the half-tolerance in phase discrepancy at MS walls
MS_WALL_PHASE_TOLERANCE = 3.14159 / 10.0

### Defines the percentage tolerance in numerics discrepancy of MS walls
### e.g. 0.01 means 1%
C_C_TOLERANCE = 0.01


#------------------------------------------------------------------------
### All following parameters concern ODE-INT in the computation
### when checking the numerics of kwall periods as compared to 
### period evolution from the W-model basepoint

### The step in the integration of periods along test paths
ODE_INT_STEP = 0.01

### Maximum tolerated absolute value of det(pf_matrix) during ode-int.
TRAJECTORY_SINGULARITY_THRESHOLD = 10**3

### Tells ode-int to stop when the integration reaches proximity of the 
### tip of a kwall, within a radius detemined by this value
ODE_INT_PROXIMITY_THRESHOLD = 0.01

### Tells ode-int to stop when the derivative of the holomorphic period
### grows larger than by this value
ODE_INT_DERIVATIVE_THRESHOLD  = 100000

#------------------------------------------------------------------------

### Defines the percentage tolerance in numerics discrepancy of 
### holomorphic periods on kwalls as compared to reference values 
### computed by PF transport from the basepoint 
###e.g. 0.01 means 1%
PERIOD_TOLERANCE = 0.01

#------------------------------------------------------------------------

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


def diagnose_kwall_network(kwn):

    print "\n\t1) Checking numerics of central charges:\n"
    for k in kwn.k_walls:
        check_c_c_numerics(k)

    print "\n\n\t2) Checking phases of central charges at MS walls:\n"
    for i in kwn.intersections:
        check_marginal_stability_condition(i)

    print "\n\n\t3) Checking positivity of DSZ pairings at intersections:\n"
    pairing_matrix = array(kwn.fibration.dsz_matrix)
    for i in kwn.intersections:
        check_dsz_positivity(i, pairing_matrix)

    # print "\n\n\t4) Checking numerics of holomorphic periods:\n"
    # for k in kwn.k_walls:
    #     check_evolved_period(kwn, k)



def check_dsz_positivity(intersection, pairing_matrix):
    ### The parents are sorted by the IntersectionPoint class
    ### in such a way that the KSWCF is K_2 K_1 = K_1 ... K_2
    ### i.e. Arg(Z_1) < Arg(Z_2)
    parents = intersection.parents
    gamma_1 = array(parents[0].charge(intersection.index_1))
    gamma_2 = array(parents[1].charge(intersection.index_2))
    m = gamma_1.dot(pairing_matrix.dot(gamma_2))

    if m < 0:
        print "NEGATIVE intersection pairing for intersection %s" % \
                                                                intersection
    else:
        print "Positive intersection pairing for intersection %s" % \
                                                                intersection




# def check_evolved_period(kwn, kwall, show_plots=False, print_values=False): 
#     """ 
#     Detemine the period of the holomorphic one form dx / y along the
#     cycle determined by 'charge'
#     Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
#     (0,1) determined at the basepoint used to compute all the monodromies.
#     Evolves such period along a path up to u2 (belonging to some kwall) and 
#     compares with the period of the kwall at that point.
#     """
    
#     fibration = kwn.fibration
#     w_model = fibration.w_model
#     length = len(kwall.coordinates)
#     charge = kwall.charge(length-1)
    
#     ### Setting up the path of integration. 
#     ### This will be an L-shaped path starting from u_0 and ending 
#     ### at the tip of the kwall in question. 
#     ### The path avoids branch cuts on the u-plane by construction.
#     u0 = w_model.base_point
#     u2 = complexify(kwall.coordinates[-1])

#     # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
#     period_data = w_model.compute_initial_periods()
#     eta_0 = period_data[0]
#     eta_prime_0 = period_data[1]
#     beta_0 = period_data[2]
#     beta_prime_0 = period_data[3]

#     def deriv(t, y):
#         u, eta, d_eta = y 

#         matrix = fibration.pf_matrix(u)
#         det_pf = abs(det(matrix))
        
#         if det_pf > TRAJECTORY_SINGULARITY_THRESHOLD:
#             singularity_check = True
#             print "\n**************\
#             \nhit the singularity threshold!\
#             \n**************\n"

#         ### A confusing point to bear in mind: here we are solving the 
#         ### ode with respect to time t, but d_eta is understood to be 
#         ### (d eta / d u), with its own  appropriate b.c. and so on!
#         u_1 = path_derivative(u, u0, u2)
#         eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
#         d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
#         return  array([u_1, eta_1, d_eta_1])

#     ode = scipy.integrate.ode(deriv)
#     ode.set_integrator("zvode")
    
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

#         if abs(u - u2) < ODE_INT_PROXIMITY_THRESHOLD or \
#                                     abs(d_eta) > ODE_INT_DERIVATIVE_THRESHOLD:
#             print "\nInterrupting PF transport of period!\n"
#             break
#         else:
#             ode.integrate(ode.t + ODE_INT_STEP)

#     u_f, eta_gamma_f, d_eta_gamma_f = ode.y 
    
#     discrepancy = abs((eta_gamma_f - k0.periods[-1]) / k0.periods[-1])

#     if show_plots == True:
#         data_plot(recorded_loci,"u along PF path")
#         data_plot(recorded_periods,"eta along PF path")
#         data_plot(recorded_d_eta,"d_eta along PF path")
#     else:
#         pass

#     if print_values == True:
#         print "u_f = %s" % u_f
#         print "eta_f = %s\n" % eta_gamma_f
#         print "to compare with:\nkwall point = %s\nkwall period = %s\n" % \
#                     (complexify(k0.coordinates[-1]), k0.periods[-1])



def check_c_c_numerics(kwall):
    #### An alternative method of computing the central charge
    central_charge_alt =  [0.0]

    for i in range(len(kwall.coordinates[:-1])):
        du = complexify(kwall.coordinates[i+1]) \
             - complexify(kwall.coordinates[i])
        eta_avg = 0.5 * (kwall.periods[i+1] + kwall.periods[i])
        c_c = complex(central_charge_alt[-1] + eta_avg * du)
        central_charge_alt.append(c_c) 

    discrepancy = abs((central_charge_alt[-1] - kwall.central_charge[-1]) / 
                                                    kwall.central_charge[-1])

    if discrepancy < C_C_TOLERANCE:
        print "Central charge numerics for kwall %s\n is OK to %s accuracy\n" \
                                                    % (kwall, C_C_TOLERANCE)



def check_marginal_stability_condition(intersection):
    kwall_1 = intersection.parents[0]
    kwall_2 = intersection.parents[1]
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    locus = intersection.locus

    Z_1 = kwall_1.central_charge[index_1]
    Z_2 = kwall_2.central_charge[index_2]

    ### DEFINE THIS NUMERICAL CONSTANT ELSEWHERE !!!
    ###
    if -1.0 * MS_WALL_PHASE_TOLERANCE \
                                    < phase(Z_1 / Z_2) \
                                    < MS_WALL_PHASE_TOLERANCE :
        # logging.debug('\nOK: the central charges of kwalls {} do align\
        #        \nat their intersection u = {}. \
        #        \nIn fact, they are:\
        #        \nZ_1 = {}\nZ_2 = {}\n'\
        #        .format([kwall_1, kwall_2], locus, Z_1, Z_2))
        print('\nOK: the central charges of kwalls {} do align\
               \nat their intersection u = {}. \
               \nIn fact, they are:\
               \nZ_1 = {}\nZ_2 = {}\n'\
               .format([kwall_1, kwall_2], locus, Z_1, Z_2))

    else:
        ### the phase discrepancy is too large to be on a MS wall
        # logging.debug('\nWARNING: the central charges of kwalls {} don-t align\
        #        \nat their intersection u = {}. \
        #        \nIn fact, they are:\
        #        \nZ_1 = {}\nZ_2 = {}\n'\
        #        .format([kwall_1, kwall_2], locus, Z_1, Z_2))
        print('\nWARNING: the central charges of kwalls {} don-t align\
               \nat their intersection u = {}. \
               \nIn fact, they are:\
               \nZ_1 = {}\nZ_2 = {}\n'\
               .format([kwall_1, kwall_2], locus, Z_1, Z_2))
    
    pass
