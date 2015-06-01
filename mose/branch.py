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
from sympy import diff, N, simplify, series, I
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
        logging.debug('\nprimary growth\n')
        self.hair = Hair(self)

        logging.debug('\nPF growth\n')
        trajectory_singularity_threshold = 10 ** 6
        ode_size_of_step = 1e-1
        
        ### SHOULD REALLY ERASE THIS PARAMETER!
        ### BUT WOULD REQUIRE TO IMPROVE THE WAY PF EVOLUTION IS HANDLED 
        ### BELOW. TO DO!
        ode_num_steps = 500000
        
        h_0 = complexify(self.hair.coordinates[0]).imag
        max_distance = 0.5 * minimum_distance(self.fibration.branch_points)
        self.hair.grow_pf(
                        # h_0=h_0,
                        # max_distance=max_distance,
                        trajectory_singularity_threshold=\
                                            trajectory_singularity_threshold,
                        ode_size_of_step=ode_size_of_step,   
                        ode_num_steps=ode_num_steps,
                        )

        d_eta_by_d_u = [(self.hair.periods[i+1] - self.hair.periods[i]) \
                    / (complexify(self.hair.coordinates[i+1]) \
                    - complexify(self.hair.coordinates[i])) \
                    for i in range(len(self.hair.coordinates)-1)]
        
        ### If the level is '10' it means it's 'debug'
        ### hence we will display this nice plot
        if logging.getLogger().getEffectiveLevel() == 10:
            data_plot(d_eta_by_d_u, "Periods along the hair")

        ###  DEIFNE THIS NUMERICAL CONSTANT ELSEWHERE !!!!
        ###
        if abs(self.hair.base_point - complexify(self.hair.coordinates[-1]))\
                                                                     > 0.01:
            logging.info('\nHair growth didn-t reach the basepoint!\n')

    def determine_positive_period(self, reference_period):
        hair_initial_period = self.hair.periods[0]
        hair_final_period = self.hair.periods[-1]
        sign = periods_relative_sign(\
                                    reference_period, \
                                    hair_final_period \
                                    )
        self.positive_period = hair_initial_period * sign

        logging.debug('\nThe hair final period is : {}\nat : {}\n'.format(\
                                hair_final_period, self.hair.coordinates[-1]))
        logging.debug('\nThe reference period is : {}\n'\
                                                    .format(reference_period))



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
        u0 = self.initial_point.locus
        u = sym.Symbol('u')
        v = sym.Symbol('v')
        x = sym.Symbol('x')

        # ### Disabling the following for this reason:
        # ### Instead of computing these for each branch point, 
        # ### we do it once and for all in the elliptic_fibration module
        # ### and then pass the result here
        # w_f = self.fibration.num_f
        # w_g = self.fibration.num_g
        # eq = sym.simplify(x ** 3 + w_f * x + w_g)
        # self.initial_point.sym_roots = sym.simplify(sym.solve(eq, x))



        # ### There is no apparent gain in speed using the approximation below
        # ### In fact, it's still a cubic equation that we need to solve,
        # ### no matter what.
        # ### But keep it here, just in case it might be relevant for numerics
        # ### with complicated fibrations...
        # # ### To speed up the computation of hair initial evolution, as well 
        # # ### as primary Kwall initial evolution, we consider a Taylor
        # # ### expansion in u, instead of the full expression.
        # # ### This will help with the computation of the roots of the equation.
        # # eq_1 = sym.simplify(eq.subs(u, u0 + v))
        # # eq_2 = series(eq_1, x=u, n=2).removeO().subs(v, u - u0)
        # # self.initial_point.sym_roots = sym.simplify(sym.solve(eq_2, x))

        e1, e2, e3 = self.fibration.sym_roots

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

            f1, f2, f3 = map(complex, map(lambda x: x.subs(u, u1), \
                                                self.fibration.sym_roots))

            roots = [f1, f2, f3]

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

        self.pf_boundary_conditions = [u1, eta_1, eta_prime_1]
    
    def grow_pf(self,
                trajectory_singularity_threshold=None,
                ode_size_of_step=None, ode_num_steps=None, **ode_kwargs):
        """ 
        Implementation of the Picard-Fuchs growth of hairs
        """

        if(ode_num_steps is None or ode_num_steps == 0):
            # Don't grow this hair, exit immediately.
            return None

        ode = scipy.integrate.ode(hair_pf_ode_f)
        ode.set_integrator("zvode", **ode_kwargs)

        y_0 = self.pf_boundary_conditions
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
            ### checking if we reached the basepoint, in that case 
            ### then ode.y would be ['stop','stop','stop']
            u, eta, d_eta = ode.y
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
    if u_1 == 'stop':
        logging.debug('Stopping hair growth')
        hair.growth_control = 'stop'
        return array([u, eta, d_eta])

    else:
        eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)    
        return  array([u_1, eta_1, d_eta_1])

        
def reference_period(branch_point): 
    """ 
    Detemine the period of the holomorphic one form dx / y along the
    cycle that pinches at a given discriminant locus (the n-th).
    Uses Picard-Fuchs evolution and initial values of the periods (1,0) and
    (0,1) determined at the basepoint used to compute all the monodromies.
    """
    
    logging.debug('\
        \n******************************************\
        \nEvolving periods for discriminant locus {}\
        \n******************************************\n'.format(\
                                                        branch_point.count))

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
    beta_0 = period_data[1]

    ### the pinching-cycle period at the starting point, 
    eta_gamma_0 = charge[0] * eta_0 + charge[1] * beta_0
    
    return eta_gamma_0


