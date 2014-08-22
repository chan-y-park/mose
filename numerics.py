def find_singularities(g2, g3):
    import sympy as sym
    u = sym.Symbol('u')
    discriminant = g2 ** 3 - 27 * g3 ** 2
    return sym.solve(discriminant, u)


def grow_primaries(singularities, g2, g3, theta, options): 
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

    primary_kwalls = []

    for u0 in singularities:
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

        for epsilon in [(+1),(-1)]:
            def eta_1(z) : 
                return epsilon*complex( 4 * (f3 - f1).subs(u, z) ** (-0.5)  * mp.ellipk( complex( N((f2 - f1) / (f3 - f1) ).subs(u,z) ) ** 0.5 ))

            eta0 = eta_1(u0)
            bc = [u0,eta0]
            t0 = 0
            y0 = array([(complex(u0)).real,(complex(u0)).imag]) 

            def deriv(y,t):
                z = y[0]+1j * y[1]
                c = complex(N(exp( 1j * ( theta + pi - phase( eta_1(z) ) ) )))
                return array([c.real, c.imag])

            time = linspace(options[0],options[1],options[2])
            y = odeint(deriv, y0, time)
            primary_kwalls.append( y )

    return primary_kwalls


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

    time = linspace(options[0],options[1],options[2])
    y = odeint(deriv, y0, time)

    return y


