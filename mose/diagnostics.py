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


def check_evolved_period(charge, d): 
    """ 
    Detemine the period of the holomorphic one form dx / y along the
    cycle determined by 'charge'
    Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
    (0,1) determined at the basepoint used to compute all the monodromies.
    Evolves such period along a path up to u2 (belonging to some kwall) and 
    compares with the period of the kwall at that point.
    """
    kwn = d['k_wall_networks'][0]
    k0 = kwn.k_walls[0]

    fibration = kwn.fibration
    w_model = fibration.w_model
    

    ### Setting up the path of integration. Since PF blows up
    ### at the discriminant locus (although the period does not)
    ### we will not get exactly to the locus.
    u0 = w_model.base_point
    # u2 = complexify(bp.hair.coordinates[-1])
    u2 = complexify(k0.coordinates[100])
    

    
    # Notation: using "eta" for the period of (1,0), using "beta" for (0,1)
    period_data = w_model.compute_initial_periods()
    eta_0 = period_data[0]
    eta_prime_0 = period_data[1]
    beta_0 = period_data[2]
    beta_prime_0 = period_data[3]

    print "\nu_0, eta_0, beta_0:\n%s\n" % [u0, eta_0, beta_0]

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
        if abs(u - u2) < 0.01 or abs(d_eta) > 100000:
            print "\nInterrupting PF transport of period!\n"
            break
        else:
            # print "time: %s" % ode.t
            ode.integrate(ode.t + dt)

    u_f, eta_gamma_f, d_eta_gamma_f = ode.y 
    print "u_f = %s" % u_f
    print "eta_f = %s\n" % eta_gamma_f

    print "to compare with:\nkwall point = %s\nkwall period = %s\n" % \
                (complexify(k0.coordinates[100]), k0.periods[100])

    data_plot(recorded_loci,"u along PF path")
    data_plot(recorded_periods,"eta along PF path")
    data_plot(recorded_d_eta,"d_eta along PF path")


    # return eta_gamma_f





######### 
# IMPLEMENT THE CHECK WITH THE FOLLOWING ROUTINE FOR CENTRAL CHARGE NUMERICS (?)
#########

        # #### An alternative method of computing the central charge
        # self.central_charge_alt =  [0.0]

        # for i in range(len(self.coordinates[:-1])):
        #     du = complexify(self.coordinates[i+1]) \
        #          - complexify(self.coordinates[i])
        #     eta_avg = 0.5 * (self.periods[i+1] + self.periods[i])
        #     c_c = complex(self.central_charge_alt[-1] + eta_avg * du)
        #     self.central_charge_alt.append(c_c) 