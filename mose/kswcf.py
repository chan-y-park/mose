import logging

from sympy.abc import x, y, t
from sympy import Subs, series, LT, degree_list, Poly, degree
from sympy.core.numbers import NaN
from numpy import array

# To enhance the speed, we will keep memory of those KSWCFs that have been 
# already computed once. 
# Entries of stored_keys are of the form [m, omega_1, omega_2]
stored_keys = []    
# Entries of stored_results are the output of KS2()
stored_results = []     

def KS2(m, omega_1, omega_2, ks_filtration_degree):    
    """
    Solves the basic KS identity: 
    K_(0,1)^{omega_2} K_(1,0)^{omega_1} = (???)
    with 
    < (1,0) , (0,1) > = m
    We use the procedure suggested in eq. (5.2) of WWC 
    (wild wall crossing paper)

    The data should be given in such a way that omega_1
    is the degeneracy of the state with smaller phase
    Arg(Z_1) < Arg(Z_2)
    """

    ### Start by checking if we already computed this
    case = position_sublist(stored_keys, [m, omega_1, omega_2])

    if not case == None:
        return stored_results[case]

    else:
        ### The reference data: the LHS of the KSWCF
        data = [ [[0,1], omega_2],   [[1,0],omega_1]]
        var_x = S(m, data, x, ks_filtration_degree)
        var_y = S(m, data, y, ks_filtration_degree)

        ### Finding the RHS
        ansatz = [ [[1,0],omega_1], [[0,1], omega_2]]

        while True:
            ans_x = S(m, ansatz, x, ks_filtration_degree)
            ans_y = S(m, ansatz, y, ks_filtration_degree)

            expr_x = t_min_degree(
                ((var_x - ans_x) / x).expand(), ks_filtration_degree
            ) 
            # print "expr_x = %s" % expr_x
            expr_y = t_min_degree(
                ((var_y - ans_y) / y).expand(), ks_filtration_degree
            ) 
            # print "expr_y = %s" % expr_y

            if expr_x == 0 and expr_y ==0:
                break

            else:
                all_monomials_x = monomials_of_degree(
                    t_degree(expr_x, ks_filtration_degree)
                )
                # print all_monomials_x
                all_coefficients_x = [Poly(expr_x).coeff_monomial(mono) 
                                        for mono in all_monomials_x]
                # print all_coefficients_x
                all_charges_x = charges_of_degree(
                    t_degree(expr_x, ks_filtration_degree)
                )
                # print all_charges_x
                all_degeneracies_x = [
                    -all_coefficients_x[i] / 
                    (m * (-1)**(m * all_charges_x[i][0] * 
                        all_charges_x[i][1]) * all_charges_x[i][1]) 
                    for i in range(len(all_coefficients_x))
                ]
                # print all_degeneracies_x

                all_monomials_y = monomials_of_degree(
                    t_degree(expr_y, ks_filtration_degree)
                )
                # print all_monomials_y
                all_coefficients_y = [Poly(expr_y).coeff_monomial(mono) 
                                        for mono in all_monomials_y]
                # print all_coefficients_y
                all_charges_y = charges_of_degree(
                    t_degree(expr_y, ks_filtration_degree)
                )
                # print all_charges_y
                all_degeneracies_y = [
                    -all_coefficients_y[i] / 
                    (- m * (-1)**(m * all_charges_y[i][0] * 
                        all_charges_y[i][1]) * all_charges_y[i][0]) 
                    for i in range(len(all_coefficients_y))
                ]
                # print all_degeneracies_y

                all_states = []
                for i in range(len(all_degeneracies_x)):
                    omega = all_degeneracies_x[i]
                    gamma = all_charges_x[i]
                    if omega == 0 or type(omega) is NaN:
                        pass
                    else:
                        all_states.append([gamma,int(omega)])
                for i in range(len(all_degeneracies_y)):
                    omega = all_degeneracies_y[i]
                    gamma = all_charges_y[i]
                    if omega == 0 or type(omega) is NaN:
                        pass
                    else:
                        all_states.append([gamma,int(omega)])

                all_states = remove_duplicates(all_states)
                # print "new states: %s" % all_states

                # Now add the newly found BPS states to the ansatz, 
                # and reorder according to charge slope
                ansatz = ansatz + all_states
                # print ansatz

                # decorate a list of states by the charge slope, 
                # adding +0.01 to avoid dividing by zero
                decorated = [
                    [(state[0][1]+0.01)/(state[0][0]+0.01), state] 
                    for i, state in enumerate(ansatz)
                ] 
                ansatz = [state for slope, state in sorted(decorated)]
                # print ansatz

        stored_keys.append([m, omega_1, omega_2])
        stored_results.append(ansatz)

        return ansatz

def progeny_2(data, dsz, ks_filtration_degree):
    #### The formatting of "data" should be as follows:
    #### [ [gamma_2 , omega_2]  ,  [gamma_1 , omega_1] ] 
    #### phase ordered from right to left, meaning that the
    #### spectrum generator should be K_2 K_1 on one side
    #### and K_1 ... K_2 on the other side.
    #### In other words, Arg(Z_1) < Arg(Z_2)

    ### turning the list into a numpy array to perform operations on it
    gamma_1 = array(data[0][0]) 
    gamma_2 = array(data[1][0])

    omega_1 = data[0][1]
    omega_2 = data[1][1]

    logging.info('computing the progeny of : %s and %s', gamma_1, gamma_2)

    pairing_matrix = array(dsz) ### turning the list into a matrix of numpy
    m = gamma_2.dot(pairing_matrix.dot(gamma_1))
    logging.info('intersection pairing : %s', m)

    if m == 0:
        spectrum = []
    elif m > 0:
        spectrum = KS2(m, omega_1, omega_2, ks_filtration_degree)
        return spectrum[1:-1]
        # the following command would return the new states in the global 
        # basis, as opposed to the parent's basis.

        # Dropping the first and last state, because they correspond to the 
        # parents.
        #return [[(state[0][0] * gamma_1 + 
        #           state[0][1] * gamma_2).tolist(), state[1]] 
        #           for state in spectrum[1:-1]] 
    elif m < 0:
        logging.info(
                """
                \n*******************************
                \nNegative intersection pairing !
                \n*******************************
                \n\n(will use the absolute value and keep going)\n
                """
            )
        spectrum = KS2(-m,omega_2,omega_1, ks_filtration_degree)
        return list(reversed(spectrum[1:-1]))
        # the following command would return the new states in the global 
        # basis, as opposed to the parent's basis

        # Dropping the first and last state, because they correspond to the 
        # parents.
        #return [[(state[0][0] * gamma_2 + 
        #           state[0][1] * gamma_1).tolist(), state[1]] 
        #           for state in spectrum[1:-1]]

def t_expand(expr, max_deg):
    """
    Takes an expression f(x,y) and computes the Taylor expansion 
    in x and y up to bi-degree fixed by max_deg.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    x + x * y
    up to bi-degree 2.
    """
    f = expr.subs([(x, t*x), (y, t*y)])
    return series(f, t, 0, max_deg).removeO().subs(t, 1)

def t_degree(expr, max_deg):
    """
    Takes an expression f(x,y) and computes the Taylor expansion 
    in x and y up to bi-degree fixed by max_deg. Then, picks the term 
    of highest bi-degree, and returns the value of the degree.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    t * x + t**2 * x * y
    up to bi-degree 2. Will return 2.
    """
    f = expr.subs([(x, t*x), (y, t*y)])
    return degree(series(f, t, 0, max_deg).removeO(), t)


def t_min_degree(expr, max_deg, large_number=1000):
    """
    Takes an expression f(x,y) and computes the Taylor expansion 
    in x and y up to bi-degree fixed by max_deg. Then, picks the term 
    of lowest bi-degree.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    t * x + t**2 * x * y
    up to bi-degree 2. Will return x.
    """
    f = expr.subs([(x, x / t), (y, y / t)])
    f_t_series = t**large_number * series(f, t, 0, max_deg).removeO() 
    leading_term = LT(f_t_series.expand(), t).subs(t, 1)
    return leading_term 


def charges_of_degree(deg):
    return [[i, deg-i] for i in range(deg+1)]


def monomials_of_degree(deg):
    return [x**i * y**(deg-i) for i in range(deg+1)]


def remove_duplicates(a):
    """
    Removes duplicates in a list, preserving the ordering of the elements.
    """
    b = list()
    for sublist in a:
        if sublist not in b:
            b.append(sublist)
    return b


def X(gamma):
    return (x ** gamma[0]) * (y ** gamma[1])


def K(m, gamma, omega, expr, ks_filtration_degree):
    """
    This implements the KS operator, in the active form.
    The charges are understood to be expressed in a basis where 
    gamma_1 = [1,0] and gamma_2 = [0,1]
    the actual pairing matrix will be specified by a single integer,
    and reads
    [[0, m], [-m, 0]] 
    """
    x_prime = x * (
        (1 + (-1)**(m*gamma[0]*gamma[1]-1) * X(gamma))**(m*gamma[1]*omega)
    )
    y_prime = y * (
        (1 + (-1)**(m*gamma[0]*gamma[1]-1) * X(gamma))**(-m*gamma[0]*omega)
    )
    return t_expand(
        expr.subs([(x, x_prime), (y, y_prime)], simultaneous=True),
        ks_filtration_degree
    )

def S(m, data, expr, ks_filtration_degree):
    """
    This implements a sequence of KS operations.
    The formatting of "data" should be as follows.
    If 
    Arg(Z_2) < Arg(Z_1)
    then give
    [[gamma_1, omega_1], [gamma_2, omega_2]] 
    and this will compute 
    K_1 K_2
    """
    temp = expr
    for i in range(len(data)):
        gamma = data[-(i+1)][0]
        omega = data[-(i+1)][1]
        temp = K(m, gamma, omega, temp, ks_filtration_degree)
    return temp


def position_sublist(lst, sublst):
    for i in range(len(lst)):
        if sublst == lst[i]:
            return i
