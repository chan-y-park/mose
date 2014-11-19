import logging
import numpy
import scipy
import sympy as sym
import cmath
#import pdb

from sympy.utilities.lambdify import lambdify
from operator import itemgetter
from cmath import exp, pi
from numpy import array, linspace
from numpy.linalg import det
from scipy.integrate import odeint
from sympy import diff, N, simplify
from sympy import mpmath as mp
from scipy import interpolate

from branch import BranchPoint
from misc import complexify, sort_by_abs, left_right, clock
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

#    count = 0

    def __init__(self, initial_charge, degeneracy, phase, parents, 
                 fibration, boundary_condition, network, color='b'):
        self.initial_charge = initial_charge
        self.degeneracy = degeneracy
        self.phase = phase
        self.parents = parents
        self.boundary_condition = boundary_condition
        self.fibration = fibration
        self.color = color
        self.network = network
        self.singular = False
        self.cuts_intersections = []
#        logging.info('Creating K-wall #%d', KWall.count)
#        KWall.count += 1

    # def __str__(self):
    #     return ('KWall info: initial charge {}, '
    #             'degeneracy {}, etc... '.format(self.initial_charge, 
    #                                             self.degeneracy))

    def get_color(self):
        return self.color
    
    def get_xcoordinates(self):
        return [z[0] for z in self.coordinates]

    def get_ycoordinates(self):
        return [z[1] for z in self.coordinates]

    def charge(self, point):
        return self.initial_charge
        # to be updated by taking into account branch-cuts, 
        # will use data in self.splittings

    def check_cuts(self):
        #print "Checking intersections with cuts, determining charges"
        self.splittings = (55, 107, 231) 
        self.local_charge = (self.initial_charge, (2,1), (0,-1), (1,1))
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
            intersections = [[br_pt, i, clock(left_right(self.coordinates,i))] \
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



    def grow_pf(self, boundary_conditions, nint_range, 
                trajectory_singularity_threshold, pf_odeint_mxstep): 
        """ 
        implementation of the growth of kwalls 
        by means of picard-fuchs ODEs 
        """
        # NOTE: the argument boundary_conditions should be passed 
        # in the following form:
        #   u_0 = boundary_conditions[0] 
        #   eta_0 = boundary_conditions[1]
        #   (d eta / du)_0 = boundary_conditions[2] 
        # They must all be given as complex numbers. 
        # (Warning: u0 must be formatted to comply)
        # To this end, when calling this function from inside evolve() of 
        # a primary K-wall, use get_primary_bc(...)

        g2 = self.fibration.g2
        g3 = self.fibration.g3
        theta = self.phase
        u = sym.Symbol('u')
        x = sym.Symbol('x')

        Delta = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
        delta = sym.simplify(3 * (g3) * diff(g2, u) - 2 * (g2) * diff(g3, u))

        
        pf_matrix = lambdify(u, 
            [
                [0, 1], 
                [
                    sym.simplify(
                        (-(3*g2*(delta**2) / (16*(Delta**2))) + 
                        ((diff(delta,u) * diff(Delta,u)) / (12*delta*Delta)) - 
                        ((diff(Delta, u, 2)) / (12*Delta)) + 
                        (((diff(Delta, u))**2) / (144 * (Delta**2) ) ) ) 
                    ), 
                    sym.simplify( 
                        ( (diff(delta,u) / delta) - (diff(Delta,u) / Delta) ) 
                    ) 
                ]
            ] 
        )

        singularity_check= False
        t0 = 0
        y0 = map(N, 
            array(
                [
                    (boundary_conditions[0]).real,
                    (boundary_conditions[0]).imag,
                    (boundary_conditions[1]).real,
                    (boundary_conditions[1]).imag,
                    (boundary_conditions[2]).real,
                    (boundary_conditions[2]).imag
                ]
            )    
        )

        # TODO: No other way but to define deriv() here?
        def deriv(y,t):
            z = y[0] + 1j * y[1]
            eta = y[2] + 1j * y[3]
            d_eta = y[4] + 1j * y[5]
            matrix = pf_matrix(z)
            det_pf = abs(det(matrix))
            # print "PF matrix = %s" % matrix
            # print "PF determinant = %s" % det_pf
            if det_pf > trajectory_singularity_threshold:
                self.singular = True

            # A confusing point to bear in mind: here we are solving the 
            # ode with respect to time t, but d_eta is understood to be 
            # (d eta / d u), with its own  appropriate b.c. and so on!
            ### NOTE THE FOLLOWING TWO OPTIONS FOR DERIVATIVE OF z
            z_1 = exp( 1j * ( theta + pi ) ) / eta
            # z_1 = exp( 1j * ( theta + pi - phase( eta ) ) )
            eta_1 = z_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
            d_eta_1 = z_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
            return array(
                    [
                        z_1.real, z_1.imag, 
                        eta_1.real, eta_1.imag, 
                        d_eta_1.real, d_eta_1.imag
                    ]
            )
            # the following rescaled matrix will not blow up at singularities, 
            # but it will not reach the singularities either..
            # return abs(det(matrix)) * array([z_1.real, z_1.imag, 
            #              eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])

        time = linspace(*nint_range)
        y = odeint(deriv, y0, time, mxstep=pf_odeint_mxstep)

        return y


class PrimaryKWall(KWall):
    """
    K-wall that starts from a branch point.
    """
    def __init__(self, initial_charge, degeneracy, phase, parents, 
                 fibration, boundary_condition, primary_nint_range,
                 network, color='b'):
        super(PrimaryKWall, self).__init__(
            initial_charge, degeneracy, phase, parents, 
            fibration, boundary_condition, color
        )
        self.initial_point = self.parents[0]
        self.network = network

        """ 
        Implementation of the ODE for evolving primary walls, 
        valid in neighborhood of an A_1 singularity. 
        """
        g2 = self.fibration.g2
        g3 = self.fibration.g3
        theta = self.phase
        u0, sign = self.boundary_condition
        u = sym.Symbol('u')
        x = sym.Symbol('x')

        eq = 4 * x ** 3 - g2 * x - g3
        e1, e2, e3 = sym.simplify(sym.solve(eq, x))
        distances = map(abs, [e1-e2, e2-e3, e3-e1])
        pair = min(enumerate(map(lambda x: x.subs(u, u0), distances)), 
                    key=itemgetter(1))[0]
        # print "PAIR %s" % pair
        # print "e1 e2 e3: %s" % map(complex, \
        #       [e1.subs(u,u0+0.01), e2.subs(u,u0+0.01), e3.subs(u,u0+0.01)])
        if pair == 0:
            # f1 = e1
            # f2 = e2
            f1, f2 = sort_by_abs(e1, e2, u0+0.01*(1+1j))
            f3 = e3
        elif pair == 1:
            # f1 = e2
            # f2 = e3
            f1, f2 = sort_by_abs(e2, e3, u0+0.01*(1+1j))
            f3 = e1
        elif pair == 2:
            # f1 = e3
            # f2 = e1
            f1, f2 = sort_by_abs(e1, e3, u0+0.01*(1+1j))
            f3 = e2

        # print "f1, f2, f3: %s " % [f1,f2,f3]

        #TODO: No other way but to define eta_1(), deriv() here?
        eta_1_part_1 = lambdify(u, 
            # simplify(sym.expand( (sign) * 4 * (f3 - f1) ** (-0.5) )),
            (sign) * 4 * (f3 - f1) ** (-0.5),
            modules="numpy"
        ) 
        eta_1_part_2 = lambdify(u, 
            # simplify(sym.expand( ((f2 - f1) / (f3 - f1)) )),
            ((f2 - f1) / (f3 - f1)),
            modules="numpy"
        ) 
        def eta_1(z):
            return complex(eta_1_part_1(z) * mp.ellipk( eta_1_part_2(z) ) / 2.0)

        #pdb.set_trace()

        # plot real & imaginary parts of eta_1
        #plot_eta(eta_1)

        t0 = 0
        y0 = array([(complex(u0)).real,(complex(u0)).imag]) 

        ### Keep this just in case. 
        ### Modifying the boundary condition for ODE int,
        ### instead of starting from u0, start from a slightly 
        ### displaced position, to avoid trouble with potential        
        ### overlapping branch cuts from other branch-points.
        # delta = 0.01
        # u1 = u0 + delta * exp(1j*(theta + pi - cmath.phase(eta_1(u0))))
        # t0 = 0
        # y0 = array([(complex(u1)).real,(complex(u1)).imag]) 

        
        def deriv(y, t):
            z = complexify(y)
            c = exp(1j*(theta + pi - cmath.phase(eta_1(z))))
            return array([c.real, c.imag])

        start, stop, num = primary_nint_range

        time = linspace(start, stop, num)

        ##########################################
        # Save results to PrimaryKWall.coordinates
        ##########################################
        self.coordinates = odeint(deriv, y0, time)

        #
        # Calculate periods at each location of coordinates
        #
        self.periods = map(
            complex, map(eta_1, map(complexify, self.coordinates))
        )

    def get_pf_boundary_condition(self):
        """ 
        Employs data of class PrimaryKWall to produce correctly 
        formatted boundary conditions for grow_pf(). 
        """
        coords = self.coordinates
        per = self.periods
        u0 = complexify(coords[-1])
        eta0 = per[-1]
        d_eta0 = ((per[-1]-per[-2]) / 
                  (complexify(coords[-1]) - complexify(coords[-2])))
        # print "boundary conditions for PF: \n %s " % [u0, eta0, d_eta0]
        return [u0, eta0, d_eta0]

    def evolve(self, nint_range, trajectory_singularity_threshold,
               pf_odeint_mxstep):
        ti, tf, nstep = nint_range
        if(nstep == 0):
            # Don't evolve, exit immediately.
            return None
        if not (isinstance(self.parents[0], BranchPoint)):
            raise TypeError('A parent of this primary K-wall '
                            'is not a BranchPoint class.')
        # For a *primary* K-wall the boundary conditions are 
        # a bit particular:
        bc = self.get_pf_boundary_condition()

        print "PF boundary conditions: %s " % bc

        pw_data_pf = self.grow_pf(bc, nint_range,
                                  trajectory_singularity_threshold,
                                  pf_odeint_mxstep)
        self.coordinates = numpy.concatenate(
            (self.coordinates, [[row[0], row[1]] for row in pw_data_pf])
        )
        self.periods = numpy.concatenate(
            (self.periods , [row[2] + 1j* row[3] for row in pw_data_pf])
        )
        self.check_cuts()


class DescendantKWall(KWall):
    """
    K-wall that starts from an intersection point.
    """
    def __init__(self, initial_charge, degeneracy, phase, parents, fibration, 
                intersection, charge_wrt_parents, color='b'):
        """
        intersection: must be an instance of the IntersecionPoint class.
        charge_wrt_parents: must be the charge relative to 
            the parents' charge basis, hence a list of length 2.
        """
        super(DescendantKWall, self).__init__(
            initial_charge, degeneracy, phase, parents, 
            fibration, 
            [], #boundary_condition 
            color
        )
        self.initial_point = intersection
        self.charge_wrt_parents = charge_wrt_parents
        self.network = parents[0].network

    def set_pf_boundary_condition(self):
        """ 
        Employs data of class DescendantKWall to produce correctly 
        formatted boundary conditions for grow_pf(). 
        """
        intersection = self.initial_point
        charge = self.charge_wrt_parents

        u_0 = intersection.locus
        parents = intersection.parents
        
        # index_1 = intersection.index_1
        # index_2 = intersection.index_2
        # path_1 = map(complexify, parents[0].coordinates)
        # path_2 = map(complexify, parents[1].coordinates)
        # periods_1 = parents[0].periods
        # periods_2 = parents[1].periods
        # eta_1 = periods_1[index_1]
        # eta_2 = periods_2[index_2]
        # d_eta_1 = ((periods_1[index_1] - periods_1[index_1-1]) / 
        #             (path_1[index_1] - path_1[index_1-1]))
        # d_eta_2 = ((periods_2[index_2] - periods_2[index_2-1]) / 
        #             (path_2[index_2] - path_2[index_2-1]))
        # # NOTE: The use of complex() is necessary here, because sometimes 
        # # the charge vector wil be deriving from an algorithm using sympy, 
        # # and will turn j's into I's...
        # eta_0 = eta_1*complex(charge[0]) + eta_2*complex(charge[1])  
        # d_eta_0 = d_eta_1*complex(charge[0]) + d_eta_2*complex(charge[1])

        # self.boundary_condition = [u_0, eta_0, d_eta_0]


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
        
        # d_eta_1 = ((periods_1[index_1] - periods_1[index_1-1]) / 
        #               (path_1[index_1] - path_1[index_1-1]))
        # d_eta_2 = ((periods_2[index_2] - periods_2[index_2-1]) / 
        #               (path_2[index_2] - path_2[index_2-1]))
        
        ### substituted the above two lines with the following, 
        ### since in some cases you get 0 / 0 = 'nan'
        
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
        self.boundary_condition = [u_0, eta_0, d_eta_0]

    def evolve(self, nint_range, trajectory_singularity_threshold,
                pf_odeint_mxstep):
        if(nstep == 0):
            # Don't evolve, exit immediately.
            return None
        if not (isinstance(self.parents[0], KWall)):
            raise TypeError('A parent of this primary K-wall '
                            'is not a KWall class.')
        self.set_pf_boundary_condition()
        pw_data_pf = self.grow_pf(self.boundary_condition, nint_range,
                                    trajectory_singularity_threshold,
                                    pf_odeint_mxstep)
        self.coordinates =  numpy.array([[row[0], row[1]]
                                        for row in pw_data_pf])
        self.periods =  numpy.array([row[2] + 1j* row[3]
                                    for row in pw_data_pf ]) 
        self.check_cuts()

def plot_eta(eta):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure()
    X = np.arange(-10, 10, 0.1)
    Y = np.arange(-10, 10, 0.1)
    X, Y = np.meshgrid(X, Y)
    Z_min = -5.0
    Z_max = 5.0
    n_levels = 20
    levels = [Z_min + (Z_max - Z_min)*i/n_levels for i in range(n_levels)]

    eta_r = np.vectorize(lambda z: complex(eta(z)).real)
    eta_i = np.vectorize(lambda z: complex(eta(z)).imag)
    Z_r = eta_r(X + 1j * Y)
    Z_i = eta_i(X + 1j * Y)

    plt.subplot(1, 2, 1, aspect=1)
    csr = plt.contourf(X, Y, Z_r, levels=levels)
    plt.colorbar(csr, shrink=.5)

    plt.subplot(1, 2, 2, aspect=1)
    csi = plt.contourf(X, Y, Z_i, levels=levels)
    plt.colorbar(csi, shrink=.5)

    plt.show()
