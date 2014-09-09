def find_singularities(g2, g3):
    """find the singularities on the Coulomb branch"""
    import sympy as sym
    from parameters import verb

    u = sym.Symbol('u')
    discriminant = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
    if verb:
        print "\ndiscriminant: %s" % discriminant
        print "\nsingularities: %s" % sym.solve(discriminant, u)
    return sym.solve(discriminant, u)

def complexify(y):
    """ complexifies an array of two reals """
    return y[0] + 1j * y[1]

def set_primary_bc(traj):
    """ employs data of class Trajectory to produce correctly formatted boundary conditions for grow_pf """
    coords = traj.coordinates
    per = traj.periods
    u0 = complexify(coords[-1])
    eta0 = per[-1]
    d_eta0 = (per[-1]-per[-2]) / (complexify(coords[-1])-complexify(coords[-2]))
    return [u0, eta0, d_eta0]

def set_bc(intersection,charge):
    """ employs data of class Trajectory to produce correctly formatted boundary conditions for grow_pf
    Arguments: (intersection, charge)
    intersection must be an instance of the IntersecionPoint class
    charge must be the charge relative to the parents' charge basis, hence a list of length 2
    """
    u0 = intersection.locus

    parents = intersection.parents
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    path_1 = map(complexify, parents[0].coordinates)
    path_2 = map(complexify, parents[1].coordinates)
    periods_1 = parents[0].periods
    periods_2 = parents[1].periods
    eta_1 = periods_1[index_1]
    eta_2 = periods_2[index_2]
    d_eta_1 = (periods_1[index_1] - periods_1[index_1-1]) / (path_1[index_1] - path_1[index_1-1])
    d_eta_2 = (periods_2[index_2] - periods_2[index_2-1]) / (path_2[index_2] - path_2[index_2-1])
    eta0 = eta_1 * charge[0] + eta_2 * charge[1]
    d_eta0 = d_eta_1 * charge[0] + d_eta_2 * charge[1]
    return [u0, eta0, d_eta0]

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
    from sympy import N
    from sympy import mpmath as mp
    from sympy import diff
    from operator import itemgetter
    from cmath import exp, pi, phase
    from numpy import array, linspace
    from scipy.integrate import odeint

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

    t0 = 0
    y0 = map(N, array([(boundary_conditions[0]).real,(boundary_conditions[0]).imag,\
        (boundary_conditions[1]).real,(boundary_conditions[1]).imag,\
        (boundary_conditions[2]).real,(boundary_conditions[2]).imag]) )

    def deriv(y,t):
        z = y[0] + 1j * y[1]
        eta = y[2] + 1j * y[3]
        d_eta = y[4] + 1j * y[5]
        matrix = pf_matrix(z)
        #### A confusing point to bear in mind: here we are solving the ode with respect to 
        #### time t, but d_eta is understood to be (d eta / d u), with its own 
        #### appropriate b.c. and so on!
        z_1 = exp( 1j * ( theta + pi - phase( eta ) ) )
        eta_1 = z_1 * (matrix[0][0] * eta + matrix[0][1] * d_eta)
        d_eta_1 = z_1 * (matrix[1][0] * eta + matrix[1][1] * d_eta)
        return array([z_1.real, z_1.imag, eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])

    time = linspace(options[0],options[1],options[2])
    y = odeint(deriv, y0, time)

    return y


#### A TEMPORARY function to find intersections, very limited and slow
def find_intersections(trajectory_1, trajectory_2):
    path_1 = trajectory_1.coordinates
    path_2 = trajectory_2.coordinates
    from numpy import linalg as LA
    from numpy import array
    distances = [[LA.norm(array(x) - array(y)) for x in path_2 ] for y in path_1 ]
    min_dist = min(map(min, distances))
    from parameters import intersection_range
    if min_dist < intersection_range:
        for i,x in enumerate(distances):
            for j,y in enumerate(x):
                if y == min_dist:
                    index_1 = i
                    index_2 = j
                    point = complexify(list( ( array(path_1[index_1]) + array(path_2[index_2]) ) / 2 ))
                    return [[point, index_1, index_2]]  
                    ##### NOTE: the structure of the answer resembles the fact that there might be more than one intersection!
                    ##### Also, it is assumed that there should be an algorithm that determines the exact intersection point
                    ##### and inserts it into the trajectories' data, as well as computes the periods there, and the indices index_1,2 should be related to this new data point
    else:
        return []


def dsz_pairing(gamma_1, gamma_2, dsz_matrix):
    import numpy as np
    return np.dot(np.dot(gamma_1,dsz_matrix),gamma_2)


