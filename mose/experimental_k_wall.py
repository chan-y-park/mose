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

from branch import BranchPoint, minimum_distance
from misc import complexify, sort_by_abs, left_right, clock, order_roots
from monodromy import charge_monodromy

class KWall(object):
    """
    The K-wall class.

    Attributes: coordinates, periods, degeneracy, phase, charge, 
    parents, boundary_condition, count (shared), singular.
    Methods: evolve, terminate (?).
    Arguments of the instance: (initial_charge, degeneracy, phase, 
    parents, boundary_condition)
    """

    def __init__(
        self, initial_charge=None, degeneracy=None, phase=None, parents=None,
        fibration=None, color='b', label=None,
        #network=None,
    ):
        self.initial_charge = initial_charge
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
        #self.network = network
        #self.label = label
        self.singular = False
        self.cuts_intersections = []

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


    def check_cuts(self):
       ## self.splittings = (55, 107, 231) 
       ## self.local_charge = (self.initial_charge, (2,1), (0,-1), (1,1))
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

       # now sort intersections according to where they happen in proper time
       # recall that the elements of cuts_intersections  are organized as
       # [..., [branch_point, t, 'ccw'] ,...]
       # where 't' is the integer of proper time at the intersection.
       self.cuts_intersections = sorted(self.cuts_intersections , \
                                       cmp = lambda k1,k2: cmp(k1[1],k2[1]))

       # print \
       # "\nK-wall %s\nintersects the following cuts at the points\n%s\n" \
       # % (self, intersections)

       # now define the lis of splitting points (for convenience) ad the 
       # list of local charges
       self.splittings = [t for br_pt, t, chi in self.cuts_intersections]
       self.local_charge = [self.initial_charge]
       for k in range(len(self.cuts_intersections)):
           branch_point = self.cuts_intersections[k][0]   # branch-point
           # t = self.cuts_intersections[k][1]       # proper time
           direction = self.cuts_intersections[k][2]     # 'ccw' or 'cw'
           charge = self.local_charge[-1]
           new_charge = charge_monodromy(charge, branch_point, direction)
           self.local_charge.append(new_charge)


    def grow_pf(self, trajectory_singularity_threshold=None,
                ode_size_of_step=None, ode_num_steps=None, **ode_kwargs):
        """ 
        implementation of the growth of kwalls 
        by means of picard-fuchs ODEs 
        """

        if(ode_num_steps is None or ode_num_steps == 0):
            # Don't grow this K-wall, exit immediately.
            return None
        
        u = sym.Symbol('u')

        g2 = self.fibration.sym_g2.subs(self.fibration.params)
        g3 = self.fibration.sym_g3.subs(self.fibration.params)
        g2_p = diff(g2, u)
        g3_p = diff(g3, u)
        g2_p_p = diff(g2_p, u)
        g3_p_p = diff(g3_p, u)

        theta = self.phase

        #
        # Now we switch from symbolic expressions to functions
        #

        g2_n = lambdify(u, g2)
        g3_n = lambdify(u, g3)
        g2_p_n = lambdify(u, g2_p)
        g3_p_n = lambdify(u, g3_p)
        g2_p_p_n = lambdify(u, g2_p_p)
        g3_p_p_n = lambdify(u, g3_p_p)
        Delta_n = lambda z: (g2_n(z) ** 3 - 27 * g3_n(z) ** 2)
        delta_n = lambda z: (3 * (g3_n(z)) * g2_p_n(z) - 2 * \
                                                (g2_n(z)) * g3_p_n(z))

        def M10(z): 
            return (\
                -18 * (g2_n(z) ** 2) * (g2_p_n(z) ** 2) * g3_p_n(z) \
                + 3 * g2_n(z) * (7 * g3_n(z) * (g2_p_n(z) ** 3) \
                + 40 * (g3_p_n(z) ** 3)) \
                + (g2_n(z) ** 3) * (-8 * g3_p_n(z) * g2_p_p_n(z) \
                + 8 * g2_p_n(z) * g3_p_p_n(z)) \
                -108 * g3_n(z) \
                * (-2 * g3_n(z) * g3_p_n(z) * g2_p_p_n(z) \
                + g2_p_n(z) \
                * ((g3_p_n(z) ** 2) + 2 * g3_n(z) * g3_p_p_n(z))) \
                ) \
                / (16 * ((g2_n(z) ** 3) -27 * (g3_n(z) ** 2)) \
                * (-3 * g3_n(z) * g2_p_n(z) + 2 * g2_n(z) * g3_p_n(z))) 
        def M11(z):
            return \
            (-3 * (g2_n(z) ** 2) * g2_p_n(z) + 54 * g3_n(z) * g3_p_n(z)) \
            / ((g2_n(z) ** 3) - (27 * g3_n(z) ** 2)) \
            + (g2_p_n(z) * g3_p_n(z) + 3 * g3_n(z) * g2_p_p_n(z) \
            - 2 * g2_n(z) * g3_p_p_n(z)) \
            / (3 * g3_n(z) * g2_p_n(z) - 2 * g2_n(z) * g3_p_n(z))

        def pf_matrix(z): 
            return [[0, 1], [M10(z), M11(z)]]

        singularity_check= False

        def deriv(t, y):
            u, eta, d_eta = y 

            matrix = pf_matrix(u)
            det_pf = abs(det(matrix))
            if det_pf > trajectory_singularity_threshold:
                self.singular = True
                self.singular_point = u

            # A confusing point to bear in mind: here we are solving the 
            # ode with respect to time t, but d_eta is understood to be 
            # (d eta / d u), with its own  appropriate b.c. and so on!
            ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF u
            u_1 = exp( 1j * ( theta + pi ) ) / eta
            # u_1 = exp( 1j * ( theta + pi - cmath.phase( eta ) ) )
            eta_1 = u_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
            d_eta_1 = u_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
            return  array([u_1, eta_1, d_eta_1])

        y_0 = self.pf_boundary_conditions
        ode = scipy.integrate.ode(deriv)
        ode.set_integrator("zvode", **ode_kwargs)
        ode.set_initial_value(y_0)

        step = 0
        i_0 = len(self.coordinates)
        if i_0 == 0:
            self.coordinates = numpy.empty((ode_num_steps, 2), dtype=float)
            self.periods = numpy.empty(ode_num_steps, complex)
        elif i_0 > 0:
            self.coordinates.resize((i_0 + ode_num_steps, 2))
            self.periods.resize(i_0 + ode_num_steps)
        while ode.successful() and step < ode_num_steps:
            u, eta, d_eta = ode.y
            self.coordinates[i_0 + step] = [u.real, u.imag]
            self.periods[i_0 + step] = eta
            ode.integrate(ode.t + ode_size_of_step)
            step += 1
        if step < ode_num_steps:
            self.coordinates.resize((i_0 + step, 2))
            self.periods.resize(i_0 + step)

        self.check_cuts()


class PrimaryKWall(KWall):
    """
    K-wall that starts from a branch point.
    """
    def __init__(
        self, initial_charge=None, degeneracy=None, phase=None, parents=None,
        fibration=None, initial_condition=None, color='k', label=None,
        #network=None,
    ):
        if not (isinstance(parents[0], BranchPoint)):
            raise TypeError('A parent of this primary K-wall '
                            'is not a BranchPoint class.')
        super(PrimaryKWall, self).__init__(
            initial_charge=initial_charge, degeneracy=degeneracy,
            phase=phase, parents=parents, fibration=fibration, color=color,
            label=label,
            #network,
        )
        self.initial_point = self.parents[0]
        #self.network = network

        """ 
        Implementation of the ODE for evolving primary walls, 
        valid in neighborhood of an A_1 singularity. 
        """
        g2 = self.fibration.num_g2
        g3 = self.fibration.num_g3
        theta = self.phase
        u0, sign = initial_condition
        # the following parameter is to get three distinct roots e_1, e_2, e_3
        delta_u0 = (10**(-5)) * (1+1j)
        u = sym.Symbol('u')
        x = sym.Symbol('x')

        eq = sym.simplify((4 * x ** 3 - g2 * x - g3).subs(u, u0 + delta_u0))
        e1, e2, e3 = sym.simplify(sym.solve(eq, x))
        distances = map(abs, [e1 - e2, e2 - e3, e3 - e1])
        pair = min(enumerate(map(lambda x: x, distances)), 
                    key=itemgetter(1))[0]
        if pair == 0:
            f1, f2 = sort_by_abs(e1, e2)
            f3 = e3
        elif pair == 1:
            f1, f2 = sort_by_abs(e2, e3)
            f3 = e1
        elif pair == 2:
            f1, f2 = sort_by_abs(e1, e3)
            f3 = e2

        f1_0 = complex(f1.subs(u, u0))
        f2_0 = complex(f2.subs(u, u0))
        f3_0 = complex(f3.subs(u, u0))

        eta_0 = (sign) * ((f3_0 - f1_0) ** (-0.5)) * pi / 2.0

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
            roots = [f1_0, f2_0, f3_0]
            segment = [u0, u1]
            try_step = order_roots(roots, segment, sign, theta)
            # check if root tracking is no longer valid
            if try_step == 0: 
                break
            else:
                [[f1_1, f2_1, f3_1], eta_1] = try_step
                eta_prime_1 = (eta_1 - eta_0) / (u1 - u0)
                u0 = u1
                eta_0 = eta_1

        last_step = i + 1
        if last_step < max_num_steps:
            self.coordinates.resize((last_step, 2))
            self.periods.resize(last_step)
            
        self.pf_boundary_conditions = [u1, eta_1, eta_prime_1]

        ### Now we need to find the correct initial charge for the K-wall.
        ### The parent branch point gives such data, but with ambiguity on the 
        ### overall sign. This is fixed by comparing the period of the 
        ### holomorphic differential along the vanishing cycle with the 
        ### corresponding period as evolved via PF from the basepoint on the 
        ### u-plane where we trivialized the charge lattice, which we used to 
        ### compute monodromies.
        parent_bp = self.parents[0]
        positive_period = parent_bp.positive_period
        positive_charge = parent_bp.charge
        kwall_period = self.periods[0]      # this period is used for evolution
        kwall_sign = ((kwall_period / positive_period).real /
                      abs((kwall_period / positive_period).real))
        kwall_charge = list(int(kwall_sign) * array(positive_charge))
        self.initial_charge = kwall_charge


class DescendantKWall(KWall):
    """
    K-wall that starts from an intersection point.
    """
    def __init__(self, initial_charge=None, degeneracy=None, phase=None,
                 parents=None, fibration=None, intersection=None, 
                 charge_wrt_parents=None, color='b', label=None):
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
            label=label
        )
        self.initial_point = intersection
        self.charge_wrt_parents = charge_wrt_parents
        #self.network = parents[0].network
        self.pf_boundary_conditions = self.get_pf_boundary_conditions()
        self.coordinates = []
        self.periods = []

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

        if intersection.index_1 > n_trials:
            index_1 = intersection.index_1
        else:
            index_1 = n_trials ## since we need to use 'index_1 - 1 ' later on
        if intersection.index_2 > n_trials:
            index_2 = intersection.index_2
        else:
            index_2 = n_trials ## since we need to use 'index_2 - 1 ' later on
        path_1 = map(complexify, parents[0].coordinates)
        path_2 = map(complexify, parents[1].coordinates)
        periods_1 = parents[0].periods
        periods_2 = parents[1].periods
        eta_1 = periods_1[index_1]
        eta_2 = periods_2[index_2]
        
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
        return [u_0, eta_0, d_eta_0]
