class KWall:
    """
    The K-wall class.

    Attributes: coordinates, periods, degeneracy, phase, charge, 
    parents, boundary_condition, count (shared), singular.
    Methods: evolve, terminate (?).
    Arguments of the instance: (initial_charge, degeneracy, phase, 
    parents, boundary_condition)
    """

    count = 0

    def __init__(self, initial_charge, degeneracy, phase, parents, 
                 boundary_condition, color='b'):
        self.degeneracy = degeneracy
        self.phase = phase
        self.parents = parents
        self.boundary_condition = boundary_condition
        self.initial_charge = initial_charge
        self.color = color
        self.initial_point = None       # to be defined below
        self.evolve()
        Trajectory.count += 1 

    def __str__(self):
        return ('Trajectory info: initial charge {}, '
                'degeneracy {}, etc... '.format(self.initial_charge, 
                                                self.degeneracy))

    def get_color(self):
        return self.color
    
    def get_xcoordinates(self):
        return [z[0] for z in self.coordinates]

    def get_ycoordinates(self):
        return [z[1] for z in self.coordinates]

    def evolve(self):
        logging.info('Evolving trajectory %d', Trajectory.count)
        if self.parents[0].__class__.__name__ == 'BranchPoint':
            # For a *primary* kwall the boundary conditions are 
            # a bit particular:
            self.initial_point = self.parents[0]
            u0, sign = self.boundary_condition
            self.coordinates, self.periods = \
                grow_primary_kwall(u0, sign, g2, g3, self.phase, 
                                    primary_options)

            # now switch to picard-fuchs evolution
            bc = set_primary_bc(self)
            pf_evolution = grow_pf(bc, g2, g3, self.phase, options)
            pw_data_pf = pf_evolution[0]
            self.singular = pf_evolution[1]
            self.coordinates = numpy.concatenate(( self.coordinates, [ [row[0], row[1]] for row in pw_data_pf ] ))
            self.periods = numpy.concatenate(( self.periods , [row[2] + 1j* row[3] for row in pw_data_pf] ))
            self.check_cuts()
        elif self.parents[0].__class__.__name__ == 'Trajectory':
            ### recall that self.boundary_conditions in this case are set by set_bc(...)
            ### and they are formatted as [u0, eta0, d_eta0, intersection]
            ### the latter being an IntersectionPoint object
            self.initial_point = self.boundary_condition[3]
            pf_evolution = grow_pf(self.boundary_condition[0:3], g2, g3, self.phase, options)
            pw_data_pf = pf_evolution[0]
            self.singular = pf_evolution[1]
            self.coordinates =  [ [row[0], row[1]] for row in pw_data_pf ]
            self.periods =  [ row[2] + 1j* row[3] for row in pw_data_pf ] 
            self.check_cuts()

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

def grow_primary_kwall(u0, sign, g2, g3, theta, primary_options): 
    """ implementation of the ODE for evolving primary walls, valid in neighborhood of an A_1 singularity """
    import sympy as sym
    from sympy import N
    from sympy import mpmath as mp
    from operator import itemgetter
    from cmath import exp, pi, phase
    from numpy import array, linspace
    from scipy.integrate import odeint

    u = sym.Symbol('u')
    x = sym.Symbol('x')

    eq = 4 * x ** 3 - g2 * x - g3
    roots = sym.solve(eq, x)
    e1 = roots[0]
    e2 = roots[1]
    e3 = roots[2]
    distances = map(abs, [e1-e2, e2-e3, e3-e1])
    def sub_u0(x): return x.subs(u,u0)

    pair = min(enumerate( map(sub_u0, distances) ), key=itemgetter(1))[0]
    if pair == 0:
        f1 = e1
        f2 = e2
        f3 = e3
    elif pair == 1:
        f1 = e2
        f2 = e3
        f3 = e1
    elif pair == 2:
        f1 = e3
        f2 = e1
        f3 = e2

    from sympy.utilities.lambdify import lambdify
    import scipy
    eta_1_part_1 = lambdify(u, sym.expand( (sign)* ( 4 * (f3 - f1) ** (-0.5) )), mp) 
    eta_1_part_2 = lambdify(u, (sym.expand( (f2 - f1) / (f3 - f1)  ) ), mp)
    def eta_1(z): 
        return eta_1_part_1(z) * mp.ellipk( eta_1_part_2(z) ) / 2

    eta0 = eta_1(u0)
    bc = [u0,eta0]
    t0 = 0
    y0 = array([(complex(u0)).real,(complex(u0)).imag]) 

    def deriv(y,t):
        z = complexify(y)
        c = exp( 1j * ( theta + pi - phase( eta_1(z) ) ) )
        return array([c.real, c.imag])

    time = linspace(primary_options[0],primary_options[1],primary_options[2])
    y = odeint(deriv, y0, time)

    return [y , map( complex, map( eta_1 , map(complexify, y) ) ) ]

def grow_pf(boundary_conditions, g2, g3, theta, options): 
    """ implementation of the growth of kwalls by means of picard-fuchs ODEs """
    import sympy as sym
    from sympy import N, diff
    from sympy import mpmath as mp
    from operator import itemgetter
    from cmath import exp, pi, phase
    from numpy import array, linspace
    from numpy.linalg import det
    from scipy.integrate import odeint
    from parameters import TRAJECTORY_SINGULARITY_THRESHOLD
    global singularity_check

    # # NOTE: the argument boundary_conditions should be passed in the following form:
    # u_0 = boundary_conditions[0] 
    # eta_0 = boundary_conditions[1]
    # (d eta / du)_0 = boundary_conditions[2] 
    # ## They must all be given as complex numbers (Warning: u0 must be formatted to comply)
    # ## to this end, when calling this function from inside evolve() of a primary Kwall, 
    # ## use set_primary_bc(...)

    u = sym.Symbol('u')
    x = sym.Symbol('x')

    Delta = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
    delta = sym.simplify(3 * (g3) * diff(g2, u) - 2 * (g2) * diff(g3, u))

    from sympy.utilities.lambdify import lambdify
    pf_matrix = lambdify(u, [ [0, 1] , [ sym.simplify(  ( -(3 * g2 * (delta ** 2) / (16 * (Delta**2)) ) + ( (diff(delta,u) * diff(Delta,u)) / (12 * delta * Delta) ) - ( (diff(Delta, u, 2)) / (12 * Delta) ) + ( ((diff(Delta, u)) ** 2) / (144 * (Delta ** 2) ) ) ) ) , sym.simplify( ( (diff(delta,u) / delta) - (diff(Delta,u) / Delta) ) ) ] ] )

    singularity_check= False
    t0 = 0
    y0 = map(N, array([(boundary_conditions[0]).real,(boundary_conditions[0]).imag,\
        (boundary_conditions[1]).real,(boundary_conditions[1]).imag,\
        (boundary_conditions[2]).real,(boundary_conditions[2]).imag]) )

    def deriv(y,t):
        global singularity_check
        z = y[0] + 1j * y[1]
        eta = y[2] + 1j * y[3]
        d_eta = y[4] + 1j * y[5]
        matrix = pf_matrix(z)
        if abs(det(matrix)) > TRAJECTORY_SINGULARITY_THRESHOLD:
            singularity_check = True

        #### A confusing point to bear in mind: here we are solving the ode with respect to 
        #### time t, but d_eta is understood to be (d eta / d u), with its own 
        #### appropriate b.c. and so on!
        z_1 = exp( 1j * ( theta + pi - phase( eta ) ) )
        eta_1 = z_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = z_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        return array([z_1.real, z_1.imag, eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])
        # the following rescaled matrix will not blow up at singularities, but it will not reach the singularities either..
        # return abs(det(matrix)) * array([z_1.real, z_1.imag, eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])

    time = linspace(options[0],options[1],options[2])
    y = odeint(deriv, y0, time, mxstep=5000000)

    return [y, singularity_check]

def set_primary_bc(traj):
    """ employs data of class Trajectory to produce correctly formatted boundary conditions for grow_pf """
    coords = traj.coordinates
    per = traj.periods
    u0 = complexify(coords[-1])
    eta0 = per[-1]
    d_eta0 = (per[-1]-per[-2]) / (complexify(coords[-1])-complexify(coords[-2]))
    return [u0, eta0, d_eta0]


