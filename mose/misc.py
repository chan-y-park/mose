import logging
import numpy
import datetime
import cmath
from scipy.optimize import fsolve
from scipy.integrate import quad as n_int
from sympy.abc import u
from sympy import mpmath as mp
# from scipy import mpmath as mp
from operator import itemgetter
from scipy.integrate import quad as n_int
from cmath import pi, exp, phase, sqrt
import matplotlib.pyplot as plt

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

    if abs(a_val) > abs(b_val):
        return a, b
    elif abs(b_val) > abs(a_val):
        return b, a
    elif abs(b_val) == abs(a_val):
        logging.info('\nCANT SORT ROOTS NEAR A DISCRIMINANT LOCUS!\n')
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
        logging.info('Cant determine direction, point doesnt belong to list!')
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
        logging.info('\nCannot read direction!\n')

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
        eta_u1 = (sign) * 4.0 * ((e3 - f1) ** (-0.5)) * \
                                    mp.ellipk( ((f2 - f1) / (e3 - f1)) )
        phase_1 = cmath.phase( \
                cmath.exp(1j * (theta + cmath.pi - cmath.phase(eta_u1))) / \
                (u1 - u0) \
                )

        ### Second possibility ###
        f1, f2 = twins[::-1]
        eta_u1 = (sign) * 4.0 * ((e3 - f1) ** (-0.5)) * \
                                    mp.ellipk( ((f2 - f1) / (e3 - f1)) )
        phase_2 = cmath.phase( \
                cmath.exp(1j * (theta + cmath.pi - cmath.phase(eta_u1))) / \
                (u1 - u0) \
                )

        if abs(phase_1) < abs(phase_2):
            e1, e2 = twins
        else:
            e1, e2 = twins[::-1]
        
        eta_u1 = (sign) * 4.0 * ((e3 - e1) ** (-0.5)) * \
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

def phase2(z):
    """
    Give the phase of a complex number z between 0 and 2pi
    """
    phase_z = phase(z)
    if phase_z < 0:
        return phase_z + 2*pi
    else:
        return phase_z

def integrand(period, *args):
    """
    2/sqrt{(x-a)(x-b)(x-c)} for x along the line (a, b) or (b, c), 
    with a choice of branch cuts.
    """
    x, a, b, c = [complex(v) for v in args]
    if period == 'A':
        # x \in (a, b)
        theta_a = phase(x-a)
        theta_b = theta_a + pi
        theta_c = phase2(x-c)
    elif period == 'B':
        # x \in (b, c)
        theta_a = phase(x-a)
        theta_b = phase(x-b)
        theta_c = theta_b + pi
    return 2/(sqrt(abs(x-a)*abs(x-b)*abs(x-c)) *
              exp(1.0j*(theta_a + theta_b + theta_c)/2.0))

# XXX: Consider a factor of 2 from (y^2 = 4x^3 + ...) vs (y^2 = x^3 + ...)
def period_A(a, b, c):
    """
    Calculate \int_b^a 2/sqrt{(x-a)(x-b)(x-c)} dx
    """
    # x = b + (a-b)*t, t \in [0, 1].
    fr = lambda t: ((a-b) * integrand('A', b+(a-b)*t, a, b, c)).real
    fi = lambda t: ((a-b) * integrand('A', b+(a-b)*t, a, b, c)).imag
    r_part, r_error = n_int(fr, 0, 1)
    i_part, i_error = n_int(fi, 0, 1)
    return (r_part + 1j*i_part, r_error, i_error)


def period_B(a, b, c):
    """
    Calculate \int_b^c 2/sqrt{(x-a)(x-b)(x-c)} dx
    """
    # x = b + (c-b)*t, t \in [0, 1].
    fr = lambda t: ((c-b) * integrand('B', b+(c-b)*t, a, b, c)).real
    fi = lambda t: ((c-b) * integrand('B', b+(c-b)*t, a, b, c)).imag
    r_part, r_error = n_int(fr, 0, 1)
    i_part, i_error = n_int(fi, 0, 1)
    return (r_part + 1j*i_part, r_error, i_error)

def periods_relative_sign(p_1, p_2):
    """
    This function determines the relative sign between 
    two numerically computed periods. 
    It also checks that they are reasonably similar, 
    it will print a message to warn if they aren't, 
    but it will let the code go on nevertheless.
    """
    trouble = ' '
    sign = 1.0

    if -1.0 * pi / 4.0 < phase(p_1 / p_2) < 1.0 * pi / 4.0:
        sign = 1.0
    elif 3.0 * pi / 4.0 < phase(p_1 / p_2) <= pi or \
       -1.0 * pi <= phase(p_1 / p_2) < -3.0 * pi / 4.0:
       sign = -1.0
    else:
        trouble += ' phase discrepancy too large, '

    if 0.5 < abs(p_1) / abs(p_2) < 2.0:
        pass
    else:
        trouble += ' modulus discrepancy too large '

    if trouble != ' ':
        logging.info('\
            \nWARNING: could not reliably determine the positive period, \
            \ndue to: {}'.format(trouble))

    return sign

def check_marginal_stability_condition(intersection):
    ### Enable the code below, once the computation of 
    ### central charges is implemented.

    kwall_1 = intersection.parents[0]
    kwall_2 = intersection.parents[1]
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    locus = intersection.locus

    Z_1 = kwall_1.central_charge[index_1]
    Z_2 = kwall_2.central_charge[index_2]

    ### DEFINE THIS NUMERICAL CONSTANT ELSEWHERE !!!
    ###
    if -1.0 * pi / 10 < phase(Z_1 / Z_2) < pi / 10 :
        logging.debug('\nOK: the central charges of kwalls {} do align\
               \nat their intersection u = {}. \
               \nIn fact, they are:\
               \nZ_1 = {}\nZ_2 = {}\n'\
               .format([kwall_1, kwall_2], locus, Z_1, Z_2))

    else:
        ### the phase discrepancy is too large to be on a MS wall
        logging.debug('\nWARNING: the central charges of kwalls {} don-t align\
               \nat their intersection u = {}. \
               \nIn fact, they are:\
               \nZ_1 = {}\nZ_2 = {}\n'\
               .format([kwall_1, kwall_2], locus, Z_1, Z_2))
    
    pass



def sort_parent_kwalls(parents, indices):
    ### Enable the code below, once the computation of 
    ### central charges is implemented.

    kwall_1 = parents[0]
    kwall_2 = parents[1]
    ### will check central charges slightly before the walls intersect
    index_1 = max(indices[0] - 10, 0)
    index_2 = max(indices[1] - 10, 0)

    Z_1 = kwall_1.central_charge[index_1]
    Z_2 = kwall_2.central_charge[index_2]

    if phase(Z_1 / Z_2) > 0:
        ### the phase of Z_1 is greater than the phase of Z_2
        return [kwall_1, kwall_2]
    else:
        ### the phase of Z_1 is smaller than the phase of Z_2
        return [kwall_2, kwall_1]

    # return parents


def data_plot(cmplx_list, title):
    """
    Plot the real and imaginary parts of 
    a list of complex numers
    """

    l = len(cmplx_list)
    r_list = [x.real for x in cmplx_list]
    i_list = [x.imag for x in cmplx_list]

    plt.figure(1)

    plt.subplot(211)
    plt.plot(r_list, "r.")
    r_delta = abs(max(r_list)-min(r_list))
    plt.axis([0.0, float(l), min(r_list) - 0.15 * r_delta, \
                             max(r_list) + 0.15 * r_delta])
    plt.ylabel("real part")
    plt.title(title)

    plt.subplot(212)
    plt.plot(i_list, "r.")
    i_delta = abs(max(i_list)-min(i_list))
    plt.axis([0.0, float(l), min(i_list) - 0.15 * i_delta, \
                             max(i_list) + 0.15 * i_delta])
    plt.ylabel("imaginary part")
    plt.title(title)
    
    plt.show()


def path_derivative(u, u_i, u_f):
    """
    This function returns the derivative du/dt 
    for the path u(t) that is to be used for 
    computing the positive period via PF transport
    FROM the basepoint.

    parameters:
    -----------

    u = point of evaluation
    u_i = the initial basepoint
    u_f = the final point (the disc. locus)

    the output:
    -----------

    it will give a posivive real derivative 
    until Re(u) reaches Re(u_f).
    Then it will give a positive imaginary
    derivative, until Im(u) reaches Im(u_f).
    In this way, an L-shaped path is produced
    """

    epsilon = 0.01

    if u_i.real < u_f.real and u_i.imag < u_f.imag:
        if u.real < u_f.real:
            ### The last factor is meant to make it slow down
            ### as it approaches the right 'height' (imaginary part)
            ### but the '+epsilon' is needed, otherwise it will go 
            ### to zero, and will not make the turn to continue evolving
            ### towards u_f.
            ### The same reasonning applies everywhere else in this function.
            return epsilon * (u_f.real - u.real)
        elif u.imag < u_f.imag:
            return epsilon * 1j * (u_f.imag - u.imag)
        else: 
            return 0

    elif u_i.real < u_f.real and u_i.imag > u_f.imag:
        if u.real < u_f.real:
            return epsilon * (u_f.real - u.real)
        elif u.imag > u_f.imag:
            return -1 * epsilon * 1j * (u.imag - u_f.imag)
        else: 
            return 0

    if u_i.real > u_f.real and u_i.imag < u_f.imag:
        if u.real > u_f.real:
            return -1 * epsilon * (u.real - u_f.real)
        elif u.imag < u_f.imag:
            return epsilon * 1j * (u_f.imag - u.imag)
        else: 
            return 0

    if u_i.real > u_f.real and u_i.imag > u_f.imag:
        if u.real > u_f.real:
            return -1 * epsilon * (u.real - u_f.real)
        elif u.imag > u_f.imag:
            return -1 * epsilon * 1j * (u.imag - u_f.imag)
        else: 
            return 0   


def path_derivative_alt(u, u_i, u_f):
    """
    This function returns the derivative du/dt 
    for the path u(t) that is to be used for 
    computing the positive period via PF transport
    TO the basepoint.
    The difference with the former function is that
    here we give priority to vertical motion, 
    and only then proceed horizontally.

    parameters:
    -----------

    u = point of evaluation
    u_i = the initial basepoint
    u_f = the final point (the disc. locus)

    the output:
    -----------

    it will give a potsivive real derivative 
    until Re(u) reaches Re(u_f).
    Then it will give a positive imaginary
    derivative, until Im(u) reaches Im(u_f).
    In this way, an L-shaped path is produced
    """

    epsilon = 0.01

    if u_i.real < u_f.real and u_i.imag < u_f.imag:
        if u.imag < u_f.imag:
            ### The last factor is meant to make it slow down
            ### as it approaches the right 'height' (imaginary part)
            ### but the '+epsilon' is needed, otherwise it will go 
            ### to zero, and will not make the turn to continue evolving
            ### towards u_f.
            ### The same reasonning applies everywhere else in this function.
            return epsilon * 1j * (u_f.imag - u.imag + epsilon)
        elif u.real < u_f.real:
            return epsilon * (u_f.real - u.real + epsilon)
        else: 
            return 'stop'

    elif u_i.real < u_f.real and u_i.imag > u_f.imag:
        if u.imag > u_f.imag:
            return epsilon * (-1j) * (u.imag - u_f.imag+ epsilon)
        elif u.real < u_f.real:
            return epsilon * (u_f.real - u.real + epsilon)
        else: 
            return 'stop'

    if u_i.real > u_f.real and u_i.imag < u_f.imag:
        if u.imag < u_f.imag:
            return epsilon * 1j * (u_f.imag - u.imag + epsilon)
        elif u.real > u_f.real:
            return -1 * epsilon * (u.real - u_f.real + epsilon)
        else: 
            return 'stop'

    if u_i.real > u_f.real and u_i.imag > u_f.imag:
        if u.imag > u_f.imag:
            return epsilon * (-1j) * (u.imag - u_f.imag + epsilon)
        elif u.real > u_f.real:
            return -1 * epsilon * (u.real - u_f.real + epsilon)
        else: 
            return 'stop'   


