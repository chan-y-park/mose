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
from sympy import diff, N, simplify, series
from sympy import mpmath as mp
from scipy import interpolate

from branch import BranchPoint, minimum_distance
from misc import complexify, sort_by_abs, left_right, clock, order_roots, \
                periods_relative_sign, data_plot, int_sign
from monodromy import charge_monodromy, flavor_charge_monodromy

### Distance from any discriminant locus 
### within which a k-wall will be deemed 'singular'
### and be cut
DISCRIMINANT_LOCI_RADIUS = 0.03

class KWall(object):
    """
    The K-wall class.

    Attributes: coordinates, periods, degeneracy, phase, charge, 
    parents, boundary_condition, singular.
    Methods: evolve, terminate (?).
    Arguments of the instance: (initial_charge, degeneracy, phase, 
    parents, boundary_condition)
    """

    def __init__(
        self, initial_charge=None, degeneracy=None, phase=None, parents=None,
        fibration=None, color='b', label=None, identifier=None,
        network=None, initial_flavor_charge=None
    ):
        self.initial_charge = initial_charge
        self.initial_flavor_charge = initial_flavor_charge
        self.degeneracy = degeneracy
        self.phase = phase
        self.parents = parents
        # NOTE: pf_boundary_conditions should be given in the following form:
        #   u_0 = boundary_conditions[0] 
        #   eta_0 = boundary_conditions[1]
        #   (d eta / du)_0 = boundary_conditions[2] 
        # They must all be given as complex numbers. 
        self.pf_boundary_conditions = None 
        self.fibration = fibration
        self.color = color
        self.network = network
        #self.label = label
        self.singular = False
        self.cuts_intersections = []
        self.splittings = None
        self.local_charge = None
        self.local_flavor_charge = None
        self.identifier = identifier
        # self.coordinates is a N-by-2 numpy array with dtype=numpy.float64.
        self.coordinate = None

    # def __str__(self):
    #     return ('KWall info: initial charge {}, '
    #             'degeneracy {}, etc... '.format(self.initial_charge, 
    #                                             self.degeneracy))

    def get_color(self):
        return self.color
    
    def get_xcoordinates(self):
        """
        Deprecated. Use get_xs()
        """
        return [z[0] for z in self.coordinates]

    def get_xs(self):
        return self.coordinates.T[0]

    def get_ycoordinates(self):
        """
        Deprecated. Use get_ys()
        """
        return [z[1] for z in self.coordinates]

    def get_ys(self):
        return self.coordinates.T[1]

    def charge(self, point):
        sp = [0]+self.splittings+[len(self.coordinates)-1]
        if point < sp[0] or point > sp[-1]:
            print "charge(pt) must be called with pt from %s to %s for this \
trajectory!" % (sp[0], sp[-1])
            return None
        else:
            for i in range(len(sp)):
                if point <= sp[i+1]:
                    return self.local_charge[i]
                    break

    def flavor_charge(self, point):
        sp = [0]+self.splittings+[len(self.coordinates)-1]
        if point < sp[0] or point > sp[-1]:
            print "charge(pt) must be called with pt from %s to %s for this \
trajectory!" % (sp[0], sp[-1])
            return None
        else:
            for i in range(len(sp)):
                if point <= sp[i+1]:
                    return self.local_flavor_charge[i]
                    break


    def check_cuts(self):
        # determine at which points the wall crosses a cut, for instance
        # (55,107,231) would mean that we change charge 3 times
        # hence self.splittings would have length, 3 while
        # self.local_charge would have length 4.
        # local charges are determined one the branch-cut data is given,
        # perhaps computed by an external function.
        disc_locus_position = [bp.locus for bp in self.fibration.branch_points]
        # the x-coordinates of the discriminant loci
        disc_x = [z.real for z in disc_locus_position]
        # parametrizing the x-coordinate of the k-wall's coordinates
        # as a function of proper time
        traj_t = numpy.array(range(len(self.coordinates)))
        traj_x = numpy.array([z[0] for z in self.coordinates])
        # traj_y = numpy.array([z[1] for z in self.coordinates])
        # f = interp1d(traj_t, traj_x, kind = 'linear')

        # all_cuts_intersections = []

        # Scan over branch cuts, see if path ever crosses one 
        # based on x-coordinates only
        for b_pt_num, x_0 in list(enumerate(disc_x)):
            g = interpolate.splrep(traj_t, traj_x - x_0, s=0)
            # now produce a list of integers corresponding to points in the 
            # k-wall's coordinate list that seem to cross branch-cuts
            # based on the x-coordinate.
            # Will get a list [i_0, i_1, ...] of intersections
            intersections = map(int, map(round, interpolate.sproot(g)))
            # removing duplicates
            intersections = list(set(intersections))
            # enforcing y-coordinate intersection criterion:
            # branch cuts extend vertically
            y_0 = self.fibration.branch_points[b_pt_num].locus.imag
            intersections = [i for i in intersections if \
                                               self.coordinates[i][1] > y_0 ]
            # adding the branch-point identifier to each intersection
            intersections = [[self.fibration.branch_points[b_pt_num], i] \
                                                   for i in intersections]
            # dropping intersections of a primary k-wall with the 
            # branch cut emanating from its parent branch-point
            # if such intersections happens at t=0
            intersections = [[br_pt, i] for br_pt, i in intersections if \
                                   not (br_pt in self.parents and i == 0)]
            # add the direction to the intersection data: either 'cw' or 'ccw'
            intersections = [[br_pt, i, clock(left_right(self.coordinates,i))]\
                           for br_pt, i in intersections]

            self.cuts_intersections += intersections
        ### Might be worth implementing an algorithm for handling 
        ### overlapping branch cuts: e.g. the one with a lower starting point 
        ### will be taken to be on the left, or a similar criterion.
        ### Still, there will be other sorts of problems, it is necessary
        ### to just rotate the u-plane and avoid such situations.

        ### now sort intersections according to where they happen in proper 
        ### time; recall that the elements of cuts_intersections are organized 
        ### as      [..., [branch_point, t, 'ccw'] ,...]
        ### where 't' is the integer of proper time at the intersection.
        self.cuts_intersections = sorted(self.cuts_intersections , \
                                       cmp = lambda k1,k2: cmp(k1[1],k2[1]))

        logging.debug(\
        '\nK-wall {}\nintersects the following cuts at the points\n{}\n'\
        .format(self.identifier, intersections))

        ### now define the lis of splitting points (for convenience) ad the 
        ### list of local charges
        self.splittings = [t for br_pt, t, chi in self.cuts_intersections]
        self.local_charge = [self.initial_charge]
        self.local_flavor_charge = [self.initial_flavor_charge]
        for k in range(len(self.cuts_intersections)):
            branch_point = self.cuts_intersections[k][0]   # branch-point
            # t = self.cuts_intersections[k][1]           # proper time
            direction = self.cuts_intersections[k][2]     # 'ccw' or 'cw'
            gauge_charge = self.local_charge[-1]
            flavor_charge = self.local_flavor_charge[-1]
            new_gauge_charge = charge_monodromy(gauge_charge, branch_point,
                                                                    direction)
            new_flavor_charge = flavor_charge_monodromy(gauge_charge, 
                                        flavor_charge, branch_point, direction)
            self.local_charge.append(new_gauge_charge)
            self.local_flavor_charge.append(new_flavor_charge)





    def grow_pf(self, trajectory_singularity_threshold=None,
                ode_size_of_step=None, ode_num_steps=None, **ode_kwargs):
        """ 
        implementation of the growth of kwalls 
        by means of picard-fuchs ODEs 
        """

        if(ode_num_steps is None or ode_num_steps == 0):
            # Don't grow this K-wall, exit immediately.
            return None

        theta = self.phase

        ode = scipy.integrate.ode(k_wall_pf_ode_f)
        ode.set_integrator("zvode", **ode_kwargs)

        y_0 = self.pf_boundary_conditions
        ### y_0 contains the following:
        ### [u_0, eta_0, d_eta_0, central_charge_0]
        ode.set_initial_value(y_0)

        matrix = self.fibration.pf_matrix
        ode.set_f_params(matrix, trajectory_singularity_threshold, theta, self)

        step = 0
        i_0 = len(self.coordinates)
        if i_0 == 0:
            self.coordinates = numpy.empty((ode_num_steps, 2),
                                           dtype=numpy.float64)
            self.periods = numpy.empty(ode_num_steps, complex)
        elif i_0 > 0:
            self.coordinates.resize((i_0 + ode_num_steps, 2))
            self.periods.resize(i_0 + ode_num_steps)
        while ode.successful() and step < ode_num_steps:
            u, eta, d_eta, c_c = ode.y
            self.coordinates[i_0 + step] = [u.real, u.imag]
            self.periods[i_0 + step] = eta
            self.central_charge.append(c_c)
            ode.integrate(ode.t + ode_size_of_step)
            step += 1
        if step < ode_num_steps:
            self.coordinates.resize((i_0 + step, 2))
            self.periods.resize(i_0 + step)

        self.check_cuts()

    def plot_periods(self):
        data_plot(self.periods, "periods of dx/y")

    def plot_central_charge(self):
        data_plot(self.central_charge, "central charges")



class PrimaryKWall(KWall):
    """
    K-wall that starts from a branch point.
    """
    def __init__(
        self, initial_charge=None, degeneracy=None, phase=None, parents=None,
        fibration=None, initial_condition=None, color='k', identifier=None,
        network=None,
    ):
        if not (isinstance(parents[0], BranchPoint)):
            raise TypeError('A parent of this primary K-wall '
                            'is not a BranchPoint class.')
        super(PrimaryKWall, self).__init__(
            initial_charge=initial_charge, degeneracy=degeneracy,
            phase=phase, parents=parents, fibration=fibration, color=color,
            identifier=identifier, network=network, initial_flavor_charge=None
        )
        self.initial_point = self.parents[0]
        self.identifier = identifier
        
        #self.network = network

        """ 
        Implementation of the ODE for evolving primary walls, 
        valid in neighborhood of an A_1 singularity. 
        """
        
        theta = self.phase
        u0, sign = initial_condition
        u = sym.Symbol('u')
        x = sym.Symbol('x')

        ### Disabling the following: roots are computed once and
        ### for all at the level of the elliptic fibration
        # w_f = self.fibration.num_f
        # w_g = self.fibration.num_g
        # eq = x ** 3 + w_f * x + w_g
        # sym_roots = sym.simplify(sym.solve(eq, x))
        
        # ### These are computed once and for all in the neighborhood of u0
        # ### when the evolution of hair needs to compute them
        # sym_roots = self.initial_point.sym_roots

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

        ellipk_period = 4.0 * ((f3_0 - f1_0) ** (-0.5)) * pi / 2.0

        ### The 'sign' variable fixes both the charge and the period 
        ### of the K-wall, relative to those of the discriminant locus 
        ### from which it emanates.
        ### The following period and charge are compatible:
        positive_period = self.parents[0].positive_period
        positive_gauge_charge = self.parents[0].gauge_charge
        positive_flavor_charge = self.parents[0].flavor_charge
        
        eta_0 = sign * positive_period
        ### The following sign will need to be used in the initial evolution
        ### of primary kwalls: it is the sign relating the actual period
        ### of the kwall (given by eta_0) to the elliptic-K-function
        ### period.
        ellipk_sign = int_sign((eta_0 / ellipk_period).real)
        self.initial_charge = list(int(round(sign)) \
                                    * array(positive_gauge_charge))
        self.initial_flavor_charge = list(int(round(sign)) \
                                    * array(positive_flavor_charge))
        

        # The initial evolution of primary kwalls is handled with an
        # automatic tuning.
        # The length of the single step is calibrated to be
        # 1/2000 th of the minimum distance between any two discriminant
        # loci.
        # The maximal number of steps is set to 400, although it may
        # be automatically truncated whenever e_1, e_2, e_3 become too
        # hard to distinguish.
        size_of_step = minimum_distance(self.fibration.branch_points)/2000.0
        max_num_steps = 400

        self.coordinates = numpy.empty((max_num_steps, 2), dtype=float)
        self.periods = numpy.empty(max_num_steps, dtype=complex) 

        for i in range(max_num_steps):
            self.coordinates[i] = [u0.real, u0.imag]
            self.periods[i] = eta_0
            if i == max_num_steps -1:
                break
            u1 = u0 + size_of_step * exp(1j*(theta + pi - cmath.phase(eta_0)))
            # u1 = u0 + size_of_step * exp(1j*(theta + pi)) /  (10 * eta_0)

            f1, f2, f3 = map(complex, map(lambda x: x.subs(u, u1), \
                                                    self.fibration.sym_roots))
            roots = [f1, f2, f3]
            # print "ROOTS %s" % roots

            segment = [u0, u1]
            try_step = order_roots(roots, segment, ellipk_sign, theta)
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
            
        ### Computing the central charge: for primary kwalls
        ### we simply integrate by hand the values of eta(u) du
        ### Later, we use Picard-Fuchs.
        
        self.central_charge = [0.0]

        for i in range(len(self.coordinates[:-1])):
            du = complexify(self.coordinates[i+1]) \
                 - complexify(self.coordinates[i])
            eta_avg = 0.5 * (self.periods[i+1] + self.periods[i])
            c_c = complex(self.central_charge[-1] + eta_avg * du)
            self.central_charge.append(c_c) 

        self.pf_boundary_conditions = [u1, eta_1, eta_prime_1, c_c]


class DescendantKWall(KWall):
    """
    K-wall that starts from an intersection point.
    """
    def __init__(self, initial_charge=None, degeneracy=None, phase=None,
                 parents=None, fibration=None, intersection=None, 
                 charge_wrt_parents=None, color='b', identifier=None,
                 network=None, initial_flavor_charge=None):
        """
        intersection: must be an instance of the IntersecionPoint class.
        charge_wrt_parents: must be the charge relative to 
            the parents' charge basis, hence a list of length 2.
        """
        if not (isinstance(parents[0], KWall)):
            raise TypeError('A parent of this primary K-wall '
                            'is not a KWall class.')
        super(DescendantKWall, self).__init__(
            initial_charge=initial_charge, degeneracy=degeneracy, 
            phase=phase, parents=parents, fibration=fibration, color=color,
            identifier=identifier, network=network, 
            initial_flavor_charge=initial_flavor_charge
        )
        self.identifier = identifier
        self.initial_point = intersection
        self.charge_wrt_parents = charge_wrt_parents
        #self.network = parents[0].network
        self.pf_boundary_conditions = self.get_pf_boundary_conditions()
        self.coordinates = []
        self.periods = []
        self.central_charge = []

    def get_pf_boundary_conditions(self):
        """ 
        Employs data of class DescendantKWall to produce correctly 
        formatted boundary conditions for grow_pf(). 
        """
        intersection = self.initial_point
        charge = self.charge_wrt_parents

        u_0 = intersection.locus
        parents = intersection.parents
       
        ### a parameter related to handling  0 / 0 cases of d_eta_1, d_eta_2
        ### specifies how far backwards one should scan in case two consecutive 
        ### steps of a kwall have same position and periods
        n_trials = 3

        if intersection.indices[0] > n_trials:
            index_1 = intersection.indices[0]
        else:
            index_1 = n_trials ## since we need to use 'index_1 - 1 ' later on
        if intersection.indices[1] > n_trials:
            index_2 = intersection.indices[1]
        else:
            index_2 = n_trials ## since we need to use 'index_2 - 1 ' later on
        path_1 = map(complexify, parents[0].coordinates)
        path_2 = map(complexify, parents[1].coordinates)
        periods_1 = parents[0].periods
        periods_2 = parents[1].periods
        eta_1 = periods_1[index_1]
        eta_2 = periods_2[index_2]
        
        central_charge_1 = parents[0].central_charge[index_1]
        central_charge_2 = parents[1].central_charge[index_2]
        
        for i in range(n_trials):
            i += 1      # start from i = 1
            if (periods_1[index_1] == periods_1[index_1-i] and 
                path_1[index_1] == path_1[index_1-i]):
                if i < n_trials: 
                    pass
                elif i == n_trials:
                    print "\n*******\nProblem: cannot compute derivative of \
                    periods for trajectory"
            else:
                d_eta_1 = ((periods_1[index_1] - periods_1[index_1-i]) / 
                            (path_1[index_1] - path_1[index_1-i]))
                break

        for i in range(n_trials):
            i += 1      # start from i = 1
            if (periods_2[index_2] == periods_2[index_2-i] 
                    and path_2[index_2] == path_2[index_2-i]):
                if i < n_trials: 
                    pass
                elif i == n_trials:
                    print "\n*******\nProblem: cannot compute derivative of \
                    periods for trajectory"
            else:
                d_eta_2 = ((periods_2[index_2] - periods_2[index_2-i]) / 
                            (path_2[index_2] - path_2[index_2-i]))
                break
        
        eta_0 = eta_1 * complex(charge[0]) + eta_2 * complex(charge[1])  
        ### The use of complex() is necessary here, because sometimes 
        # the charge vector wil be deriving from an algorithm using sympy, 
        # and will turn j's into I's...
        d_eta_0 = d_eta_1 * complex(charge[0]) + d_eta_2 * complex(charge[1])  
        ### The use of complex() is necessary here, because sometimes 
        # the charge vector wil be deriving from an algorithm using sympy, 
        # and will turn j's into I's...

        ### Now we determine the initial central charge for the new kwall
        c_c_0 = central_charge_1 * complex(charge[0]) \
              + central_charge_2 * complex(charge[1])  

        return [u_0, eta_0, d_eta_0, c_c_0]


def k_wall_pf_ode_f(t, y, pf_matrix, trajectory_singularity_threshold, theta, \
                    kwall):
    u, eta, d_eta, c_c = y 
    matrix = pf_matrix(u)

    det_pf = abs(det(matrix))
    if kwall.__class__.__name__ == 'PrimaryKWall':
        parent_bp = kwall.parents[0]
        disc_loci = [bp.locus for bp in kwall.network.fibration.branch_points if bp!=parent_bp]
    else:
        disc_loci = [bp.locus for bp in kwall.network.fibration.branch_points]
    minimum_disc_loc_distance = min([abs(u - x) for x in disc_loci])
    
    if kwall.singular == False:
        ### The above is necessary, otherwise it will keep evolving the kwall
        if det_pf > trajectory_singularity_threshold \
                or minimum_disc_loc_distance < DISCRIMINANT_LOCI_RADIUS:
            # if minimum_disc_loc_distance < DISCRIMINANT_LOCI_RADIUS:
            #     print 'kwal too close to discrminant locus at u = %s' % u
            kwall.singular = True
            kwall.singular_point = u

    if det_pf > trajectory_singularity_threshold \
                or minimum_disc_loc_distance < DISCRIMINANT_LOCI_RADIUS:
        # A confusing point to bear in mind: here we are solving the 
        # ode with respect to time t, but d_eta is understood to be 
        # (d eta / d u), with its own  appropriate b.c. and so on!
        ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF u
        u_1 = exp( 1j * ( theta + pi ) ) / eta
        # u_1 = exp( 1j * ( theta + pi - cmath.phase( eta ) ) )
        eta_1 = 0.0
        d_eta_1 = 0.0
        c_c_1 = u_1 * eta
        return  array([u_1, eta_1, d_eta_1, c_c_1])

    else:
        # A confusing point to bear in mind: here we are solving the 
        # ode with respect to time t, but d_eta is understood to be 
        # (d eta / d u), with its own  appropriate b.c. and so on!
        ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF u
        u_1 = exp( 1j * ( theta + pi ) ) / eta
        # u_1 = exp( 1j * ( theta + pi - cmath.phase( eta ) ) )
        eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        c_c_1 = u_1 * eta
        return  array([u_1, eta_1, d_eta_1, c_c_1])


