import logging
import numpy
import datetime
import cmath
from scipy.optimize import fsolve
from sympy.abc import u
from sympy import mpmath as mp
from operator import itemgetter
from scipy.integrate import quad as n_int
from cmath import pi, exp, phase, sqrt


def complexify(y):
    """ complexifies an array of two reals """
    return y[0] + 1j * y[1]

def dsz_pairing(gamma_1, gamma_2, dsz_matrix):
    return numpy.dot(numpy.dot(gamma_1, dsz_matrix), gamma_2)

def formatted_date_time():
    today = datetime.date.today()
    now = datetime.datetime.now().time().strftime("%H.%M")
    return str(today) + '-' + str(now)

def sort_by_abs(a, b, u0):
    a_val = complex(a.subs(u, u0))
    b_val = complex(b.subs(u, u0))
    # print a, b
    # print a_val, b_val

    if abs(a_val) > abs(b_val):
        return a, b
    elif abs(b_val) > abs(a_val):
        return b, a
    elif abs(b_val) == abs(a_val):
        print "\nCANT SORT ROOTS NEAR A DISCRIMINANT LOCUS!\n"
        return a, b

def left_right(list, point):
    """
    given the list 
    [..., [x, y], ...]
    and a point in the list (specified by the corresponding integer),
    determines whether x increases or decreases at that point, 
    returning repsectively 'left' or 'right'
    """
    if point > len(list)-1:
        print "Can't determine direction, point doesn't belong to list!"
    elif point > 0:
        if list[point-1][0] < list[point][0]:
            return 'right'
        else:
            return 'left'
    elif point == 0:
        if list[point][0] < list[point+1][0]:
            return 'right'
        else:
            return 'left'

def clock(direction):
    if direction == 'left':
        return 'ccw'
    elif direction == 'right':
        return 'cw'
    else:
        print "\nCannot read direction!\n"

def is_list(p): 
    return isinstance(p, list)

def deep_reverse(mylist):
    result = []
    for e in mylist:
        if isinstance(e, list):
            result.append(deep_reverse(e))
        else:
            result.append(e)
    result.reverse()
    return result

def order_roots(roots, segment, sign, theta):
    """
    Determine the roots along a trajectory segment.
    This is done by comparing the direction of the segment
    in the complex plane with that determined by the 
    evolution equation using the roots in either order.
    The right one is chosen by similarity with the direction 
    of the segment.
    e1 and e2 are the colliding roots, e3 is the far one,
    they are all expected to be evaluated at u1.
    """
    e1_0, e2_0, e3_0 = roots
    u0, u1 = segment
    distances = map(abs, [e1_0 - e2_0, e2_0 - e3_0, e3_0 - e1_0])

    # the pair e1, e2
    pair = min(enumerate(distances), key=itemgetter(1))[0]
    if pair == 0:
        twins = [e1_0, e2_0]
        e3 = e3_0
    elif pair == 1:
        twins = [e3_0, e2_0]
        e3 = e1_0
    elif pair == 2:
        twins = [e3_0, e1_0]
        e3 = e2_0

    if abs(twins[0] - twins[1]) / abs(twins[0] - e3) < 0.3:
        ### First possibility ###
        f1, f2 = twins
        eta_u1 = (sign) * 4 * ((e3 - f1) ** (-0.5)) * \
                                    mp.ellipk( ((f2 - f1) / (e3 - f1)) )
        phase_1 = cmath.phase( \
                    cmath.exp(1j * (theta + cmath.pi - cmath.phase(eta_u1))) / \
                    (u1 - u0) \
                    )

        ### Second possibility ###
        f1, f2 = twins[::-1]
        eta_u1 = (sign) * 4 * ((e3 - f1) ** (-0.5)) * \
                                    mp.ellipk( ((f2 - f1) / (e3 - f1)) )
        phase_2 = cmath.phase( \
                    cmath.exp(1j * (theta + cmath.pi - cmath.phase(eta_u1))) / \
                    (u1 - u0) \
                    )

        if abs(phase_1) < abs(phase_2):
            e1, e2 = twins
        else:
            e1, e2 = twins[::-1]
        
        eta_u1 = (sign) * 4 * ((e3 - e1) ** (-0.5)) * \
                                    mp.ellipk( ((e2 - e1) / (e3 - e1)) )

        return [[e1, e2, e3], complex(eta_u1)]
    else:
        return 0


def cut_singular_k_wall(k_wall):
    periods = k_wall.periods
    coordinates = k_wall.coordinates
    i_0 = 0
    epsilon_0 = 10
    for i, z in enumerate(coordinates):
        epsilon = abs(complexify(z) - k_wall.singular_point)
        if epsilon < epsilon_0:
            epsilon_0 = epsilon
            i_0 = i
        else:
            break

    k_wall.coordinates = coordinates[0:i_0]
    k_wall.periods = periods[0:i_0]


def integrand_A(*args):
    x, a, b, c = [complex(v) for v in args]
    arg_c = (phase(x-c)+pi)/2.0
    return -1/(sqrt(abs(x-a)*abs(x-b)*abs(x-c))*(1j)*exp(1j*arg_c))


def integrand_B(*args):
    x, a, b, c = [complex(v) for v in args]
    arg_a = phase((x-a)/(a-b))/2
    arg_b = phase((x-b)/(a-b))/2
    arg_c = (phase(x-c)+pi)/2.0
    return 1/(sqrt(abs(x-a)*abs(x-b)*abs(x-c)) * 
              exp(1.0j*(arg_a + arg_b + arg_c)))


def period_A(a, b, c):
    fr = lambda t: ((a-b) * integrand_A(b+(a-b)*t, a, b, c)).real
    fi = lambda t: ((a-b) * integrand_A(b+(a-b)*t, a, b, c)).imag
    r_part, r_error = n_int(fr, 0, 1)
    i_part, i_error = n_int(fi, 0, 1)
    return (r_part + 1j*i_part, r_error, i_error)


def period_B(a, b, c):
    fr = lambda t: ((c-b) * integrand_B(b+(c-b)*t, a, b, c)).real
    fi = lambda t: ((c-b) * integrand_B(b+(c-b)*t, a, b, c)).imag
    r_part, r_error = n_int(fr, 0, 1)
    i_part, i_error = n_int(fi, 0, 1)
    return (r_part + 1j*i_part, r_error, i_error)
