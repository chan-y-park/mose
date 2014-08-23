def find_singularities(g2, g3):
    import sympy as sym
    u = sym.Symbol('u')
    discriminant = g2 ** 3 - 27 * g3 ** 2
    return sym.solve(discriminant, u)

def complexify(y):
    return y[0] + 1j * y[1]

def set_bc(traj):
    coords = traj.coordinates
    per = traj.periods
    u0 = complexify(coords[-1])
    eta0 = per[-1]
    d_eta0 = (per[-1]-per[-2]) / (complexify(coords[-1])-complexify(coords[-2]))
    return [u0, eta0, d_eta0]

def grow_primary_kwall(u0, sign, g2, g3, theta, options): 
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

    def eta_1(z) : 
        return (sign)*complex( 4 * (f3 - f1).subs(u, z) ** (-0.5)  * mp.ellipk( complex( N((f2 - f1) / (f3 - f1) ).subs(u,z) ) ** 0.5 ))

    eta0 = eta_1(u0)
    bc = [u0,eta0]
    t0 = 0
    y0 = array([(complex(u0)).real,(complex(u0)).imag]) 

    def deriv(y,t):
        z = y[0]+1j * y[1]
        c = complex(N(exp( 1j * ( theta + pi - phase( eta_1(z) ) ) )))
        return array([c.real, c.imag])

    new_options = [0.0, 0.5,10]
    time = linspace(new_options[0],new_options[1],new_options[2])

    #time = linspace(options[0],options[1],options[2])
    y = odeint(deriv, y0, time)

    return [y ,map(complex, map(N, map( eta_1 , map(complexify, y) ) )  ) ]


def grow_pf(boundary_conditions, g2, g3, theta, options): 
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
    # ## use set_bc(...)

    u = sym.Symbol('u')
    x = sym.Symbol('x')

    # eq = 4 * x ** 3 - g2 * x - g3
    # roots = sym.solve(eq, x)
    # e1 = roots[0]
    # e2 = roots[1]
    # e3 = roots[2]
    # distances = map(abs, [e1-e2, e2-e3, e3-e1])
    # def sub_u0(x): return x.subs(u,u0)

    # pair = min(enumerate( map(sub_u0, distances) ), key=itemgetter(1))[0]
    # if pair == 0:
    #     f1 = e1
    #     f2 = e2
    #     f3 = e3
    # elif pair == 1:
    #     f1 = e2
    #     f2 = e3
    #     f3 = e1
    # elif pair == 2:
    #     f1 = e3
    #     f2 = e1
    #     f3 = e2

    Delta = g2 ** 3 - 27 * g3 ** 2
    delta = 3 * g3 * diff(g2, u) - 2 * g2 * diff(g3, u)

    def pf_matrix(z): 
        return [ [0, 1] , [ complex( ( -(3 * g2 * delta**2) / (16 * Delta**2) + (diff(delta,u) * diff(Delta,u)) / (12 * delta * Delta) - (diff(Delta, u, 2)) / (12 * Delta) + (diff(Delta, u) ** 2) / (144 * Delta ** 2 ) ).subs(u, z) ) , complex( (diff(delta) / delta - diff(Delta) / Delta).subs(u, z) ) ] ]

    #eta0 = eta_1(u0)
    #bc = [u0,eta0]
    t0 = 0
    y0 = array([(boundary_conditions[0]).real,(boundary_conditions[0]).imag,\
        (boundary_conditions[1]).real,(boundary_conditions[1]).imag,\
        (boundary_conditions[2]).real,(boundary_conditions[2]).imag]) 

    def deriv(y,t):
        z = y[0] + 1j * y[1]
        eta = y[2] + 1j * y[3]
        d_eta = y[4] + 1j * y[5]
        z_1 = complex(N(exp( 1j * ( theta + pi - phase( eta ) ) )))
        eta_1 = pf_matrix(z)[0][0] * eta + pf_matrix(z)[0][1] * d_eta
        d_eta_1 = pf_matrix(z)[1][0] * eta + pf_matrix(z)[1][1] * d_eta
        return array([z_1.real, z_1.imag, eta_1.real, eta_1.imag, d_eta_1.real, d_eta_1.imag])

    time = linspace(options[0],options[1],options[2])
    y = odeint(deriv, y0, time)

    return y

