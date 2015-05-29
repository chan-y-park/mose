import logging
import numpy
import scipy
import sympy as sym
import cmath
import pdb

from sympy.utilities.lambdify import lambdify
from operator import itemgetter
from cmath import exp, pi
from numpy import array, linspace
from numpy.linalg import det
from scipy.integrate import odeint
from sympy import diff, N, simplify
from sympy import mpmath as mp
from scipy import interpolate

from misc import complexify, sort_by_abs, left_right, clock, order_roots, \
                periods_relative_sign, data_plot, path_derivative_alt
from monodromy import charge_monodromy

class BranchPoint:
    """The BranchPoint class.

    Attributes: locus, charge
    Arguments: locus, charge
    """

    count = 0

    def __init__(
        self, 
        locus=None,
        positive_charge=None,
        monodromy_matrix=None,
        identifier=None,
        fibration=None
    ):
        self.charge = positive_charge
        self.positive_period = None
        self.locus = locus
        self.count = BranchPoint.count
        self.fibration = fibration
        ### Old version: this would use the object itself 
        ### as the ancestor in the genealogy tree.
        ### It's a nice idea, but has some tension with
        ### multiprocessing.
        # self.genealogy = self
        ### New version: use a string as identifier.
        ### This should be more robust against issues
        ### arising from multiprocessing.
        self.genealogy = identifier
        
        self.monodromy_matrix = monodromy_matrix
        BranchPoint.count += 1

    def __str__(self):
        return 'Branch point info: charge %s, locus %s ' % \
            (self.charge, self.locus)

    def grow_hair(self):
        # print "\nprimary growth\n" 
        self.hair = Hair(self)

        # print "\nPF growth\n" 
        trajectory_singularity_threshold = 10 ** 6
        ode_size_of_step = 1e-1
        
        ### SHOULD REALLY ERASE THIS PARAMETER!
        ### BUT WOULD REQUIRE TO IMPROVE THE WAY PF EVOLUTION IS HANDLED 
        ### BELOW. TO DO!
        ode_num_steps = 100000
        
        h_0 = complexify(self.hair.coordinates[0]).imag
        max_distance = 0.5 * minimum_distance(self.fibration.branch_points)
        self.hair.grow_pf(
                        # h_0=h_0,
                        # max_distance=max_distance,
                        trajectory_singularity_threshold=trajectory_singularity_threshold,
                        ode_size_of_step=ode_size_of_step,   
                        ode_num_steps=ode_num_steps,
                        )

        d_eta_by_d_u = [(self.hair.periods[i+1] - self.hair.periods[i]) \
                    / (complexify(self.hair.coordinates[i+1]) \
                    - complexify(self.hair.coordinates[i])) \
                    for i in range(len(self.hair.coordinates)-1)]
        data_plot(d_eta_by_d_u, "Periods along the hair")

        if abs(self.hair.base_point - complexify(self.hair.coordinates[-1])) > 0.01:
            print "\nHair growth didn't reach the basepoint!\n"

    def determine_positive_period(self, reference_period):
        hair_initial_period = self.hair.periods[0]
        hair_final_period = self.hair.periods[-1]
        sign = periods_relative_sign(\
                                    reference_period, \
                                    hair_final_period \
                                    )
        self.positive_period = hair_initial_period * sign

        print "\nThe hair final period is : %s\nat : %s\n" % \
                                (hair_final_period, self.hair.coordinates[-1])
        print "\nThe reference period is : %s\n" % reference_period



# class BranchCut:
#     """
#     The BranchCut class.
#     Attributes: locus, charge, monodromy_matrix
#     """
#     count = 0

#     def __init__(self, branch_point):
#         self.charge = branch_point.charge
#         self.monodromy_matrix = branch_point.monodromy_matrix
#         self.locus = branch_point.locus
#         BranchCut.count += 1


def minimum_distance(branch_points):
    loci = [bpt.locus for bpt in branch_points]
    min_dist = abs(loci[0]-loci[1])
    for i, z_i in enumerate(loci):
        for j, z_j in list(enumerate(loci))[i+1:]:
            dist = abs(z_i - z_j)
            if dist < min_dist:
                min_dist = dist
    return min_dist


class Hair:
    """
    Hair that starts from a branch point, and grow downwards.
    They carry periods of the holomorphic 1-form, somewhat 
    similar to the PrimaryKWall class.
    """
    def __init__(self, parent):
        self.initial_point = parent
        self.fibration = parent.fibration
        self.base_point = self.fibration.w_model.base_point
        self.growth_control = 'fine'

        """ 
        Implementation of the ODE for evolving the hair, 
        valid in neighborhood of an A_1 singularity. 
        """
        w_f = self.fibration.num_f
        w_g = self.fibration.num_g

        u0 = self.initial_point.locus
        # print "LOCUS = %s " % u0
        u = sym.Symbol('u')
        x = sym.Symbol('x')

        eq = sym.simplify(x ** 3 + w_f * x + w_g)
        
        sym_roots = sym.simplify(sym.solve(eq, x))
        e1, e2, e3 = sym_roots

        distances = map(abs, [e1-e2, e2-e3, e3-e1])
        pair = min(enumerate(map(lambda x: x.subs(u, u0), distances)), 
                    key=itemgetter(1))[0]
        if pair == 0:
            f1, f2 = sort_by_abs(e1, e2, u0+(10**(-5)) * (1+1j))
            f3 = e3
        elif pair == 1:
            f1, f2 = sort_by_abs(e2, e3, u0+(10**(-5)) * (1+1j))
            f3 = e1
        elif pair == 2:
            f1, f2 = sort_by_abs(e1, e3, u0+(10**(-5)) * (1+1j))
            f3 = e2

        f1_0 = complex(f1.subs(u, u0))
        f2_0 = complex(f2.subs(u, u0))
        f3_0 = complex(f3.subs(u, u0))

        eta_0 = 4.0 * ((f3_0 - f1_0) ** (-0.5)) * pi / 2.0

        ### The initial evolution of hairs is handled with an automatic tuning.
        ### The length of the single step is calibrated to be 1/2000 th of 
        ### the minimum distance between any two discriminant loci.
        ### The maximal number of steps is set to 400, although it may
        ### be automatically truncated whenever e_1, e_2, e_3 become too
        ### hard to distinguish.

        ### !!!!!
        ### NOTE THESE NUMERICAL CONSTANTS, THEY SHOULD BE THE SAME AS FOR 
        ### PRIMARY KWALLS, SO WHY DON'T WE DEFINE THEM ELSEWHERE?
        ### !!!!!
        size_of_step = minimum_distance(self.fibration.branch_points)/5000.0
        max_num_steps = 400

        self.coordinates = numpy.empty((max_num_steps, 2), dtype=float)
        self.periods = numpy.empty(max_num_steps, dtype=complex) 

        for i in range(max_num_steps):
            self.coordinates[i] = [u0.real, u0.imag]
            self.periods[i] = eta_0
            if i == max_num_steps -1:
                break
            
            ### hairs grow downwards
            u1 = u0 + size_of_step * (-1j)

            f1, f2, f3 = map(lambda x: x.subs(u, u1), sym_roots)
            roots = [complex(f1), complex(f2), complex(f3)]
            # print "ROOTS %s" % roots

            segment = [u0, u1]
            ### setting period_sign = +1 and theta such that the
            ### growth points downwards
            theta = -pi + cmath.phase(eta_0) - pi / 2
            try_step = order_roots(roots, segment, 1, theta)
            # check if root tracking is no longer valid
            if try_step == 0: 
                break
            else:
                [roots, eta_1] = try_step
                eta_prime_1 = (eta_1 - eta_0) / (u1 - u0)
                u0 = u1
                eta_0 = eta_1

        last_step = i + 1
        if last_step < max_num_steps:
            self.coordinates.resize((last_step, 2))
            self.periods.resize(last_step)
            
        # print "HAIR PERIODS\n%s" % self.periods
        # print "HAIR COORDINATES\n%s" % self.coordinates

        self.pf_boundary_conditions = [u1, eta_1, eta_prime_1]
    
    #### This is the PF evolution used when we 
    #### grow both from the basepoint and from each of the disc loci
    ####
#     def grow_pf(self, h_0=None, max_distance=None, 
#                 trajectory_singularity_threshold=None,
#                 ode_size_of_step=None, ode_num_steps=None, **ode_kwargs):
#         """ 
#         Implementation of the Picard-Fuchs growth of hairs
#         """

#         if(ode_num_steps is None or ode_num_steps == 0):
#             # Don't grow this K-wall, exit immediately.
#             return None

#         ode = scipy.integrate.ode(hair_pf_ode_f)
#         ode.set_integrator("zvode", **ode_kwargs)

#         y_0 = self.pf_boundary_conditions
#         # print "Boundary conditions for PF evolution: \n%s" % y_0
#         ### y_0 contains the following:
#         ### [u_0, eta_0, d_eta_0]
#         ode.set_initial_value(y_0)

#         matrix = self.fibration.pf_matrix

#         step = 0
#         i_0 = len(self.coordinates)
#         u = y_0[0]
#         if i_0 == 0:
#             self.coordinates = numpy.empty((ode_num_steps, 2), dtype=float)
#             self.periods = numpy.empty(ode_num_steps, complex)
#         elif i_0 > 0:
#             self.coordinates.resize((i_0 + ode_num_steps, 2))
#             self.periods.resize(i_0 + ode_num_steps)
#         while ode.successful() and step < ode_num_steps \
#                                         and ((h_0 - u.imag) < max_distance):
#             u, eta, d_eta = ode.y
#             ode.set_f_params(matrix, trajectory_singularity_threshold, self)
#             self.coordinates[i_0 + step] = [u.real, u.imag]
#             self.periods[i_0 + step] = eta
#             ode.integrate(ode.t + ode_size_of_step)
#             step += 1
#         if step < ode_num_steps:
#             self.coordinates.resize((i_0 + step, 2))
#             self.periods.resize(i_0 + step)

# def hair_pf_ode_f(t, y, pf_matrix, trajectory_singularity_threshold, hair):
#     u, eta, d_eta = y 
#     matrix = pf_matrix(u)

#     det_pf = abs(det(matrix))
#     if det_pf > trajectory_singularity_threshold:
#         hair.singular = True
#         hair.singular_point = u

#     # A confusing point to bear in mind: here we are solving the 
#     # ode with respect to time t, but d_eta is understood to be 
#     # (d eta / d u), with its own  appropriate b.c. and so on!
#     ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF u
#     u_1 = exp( 1j * ( -pi / 2 ) ) / abs(eta)
#     # u_1 = exp( 1j * ( theta + pi ) ) / eta
#     # u_1 = exp( 1j * ( theta + pi - cmath.phase( eta ) ) )
#     eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
#     d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
#     return  array([u_1, eta_1, d_eta_1])


    #### This is the PF evolution used when we 
    #### grow only from each of the disc loci
    #### but not from the basepoint
    ####
    def grow_pf(self,
                trajectory_singularity_threshold=None,
                ode_size_of_step=None, ode_num_steps=None, **ode_kwargs):
        """ 
        Implementation of the Picard-Fuchs growth of hairs
        """

        if(ode_num_steps is None or ode_num_steps == 0):
            # Don't grow this K-wall, exit immediately.
            return None

        ode = scipy.integrate.ode(hair_pf_ode_f)
        ode.set_integrator("zvode", **ode_kwargs)

        y_0 = self.pf_boundary_conditions
        # print "Boundary conditions for PF evolution: \n%s" % y_0
        ### y_0 contains the following:
        ### [u_0, eta_0, d_eta_0]
        ode.set_initial_value(y_0)

        matrix = self.fibration.pf_matrix

        step = 0
        i_0 = len(self.coordinates)
        u = y_0[0]
        if i_0 == 0:
            self.coordinates = numpy.empty((ode_num_steps, 2), dtype=float)
            self.periods = numpy.empty(ode_num_steps, complex)
        elif i_0 > 0:
            self.coordinates.resize((i_0 + ode_num_steps, 2))
            self.periods.resize(i_0 + ode_num_steps)
        while ode.successful() and u.real > self.base_point.real \
                and step < ode_num_steps and self.growth_control != 'stop':
            u, eta, d_eta = ode.y
            ### checking if we reached the basepoint, in that case 
            ### then ode.y would be ['stop','stop','stop']
            ode.set_f_params(matrix, trajectory_singularity_threshold, self)
            self.coordinates[i_0 + step] = [u.real, u.imag]
            self.periods[i_0 + step] = eta
            ode.integrate(ode.t + ode_size_of_step)
            step += 1
        if step < ode_num_steps:
            self.coordinates.resize((i_0 + step, 2))
            self.periods.resize(i_0 + step)


    
def hair_pf_ode_f(t, y, pf_matrix, trajectory_singularity_threshold, hair):
    u_i = complexify(hair.coordinates[0])
    u_f = hair.base_point
    u, eta, d_eta = y 
    matrix = pf_matrix(u)

    det_pf = abs(det(matrix))
    if det_pf > trajectory_singularity_threshold:
        hair.singular = True
        hair.singular_point = u

    # A confusing point to bear in mind: here we are solving the 
    # ode with respect to time t, but d_eta is understood to be 
    # (d eta / d u), with its own  appropriate b.c. and so on!
    ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF u
    u_1 = path_derivative_alt(u, u_i, u_f)
    eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
    d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
    
    ### when the derivative is zero, it's time to stop!
    if u_1 == 0:
        print "Stopping hair growth"
        hair.growth_control = 'stop'

    return  array([u_1, eta_1, d_eta_1])

        

