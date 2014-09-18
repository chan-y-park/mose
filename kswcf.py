from sympy.abc import x, y, t
from sympy import Subs, series, LT, degree_list, Poly, degree
from sympy.core.numbers import NaN
from parameters import ks_filtration_degree as max_deg
from parameters import dsz_matrix as dsz
from parameters import verb
from numpy import array

#### To enhance the speed, we will keep memory of those KSWCFs that have been 
#### already computed once. 
stored_keys = []    ## Its entries are of the form [m, omega_1, omega_2]
stored_results = []     ## Its entries are the output of KS2(m, omega_1, omega_2)

def progeny_2(data):
    #### The formatting of "data" should be as follows:
    #### [ [gamma_1 , omega_1]  ,  [gamma_2 , omega_2] ] 
    #### phase ordered from right to left
    gamma_1 = array(data[0][0]) ### turning the list into a numpy array to perform operations on it
    gamma_2 = array(data[1][0]) ### same as above
    omega_1 = data[0][1]
    omega_2 = data[1][1]

    if verb:
        print "\ncomputing the progeny of : %s and %s" % (gamma_1, gamma_2)

    pairing_matrix = array(dsz) ### turning the list into a matrix of numpy
    m = gamma_1.dot(pairing_matrix.dot(gamma_2))
    if verb:
        print "\nintersection pairing : %s" % m

    if m == 0:
        spectrum = []
    elif m > 0:
        spectrum = KS2(m,omega_1,omega_2)
        return spectrum[1:-1]
        ### the following command would return the new states in the global basis, as opposed to the parent's basis
        # return [[(state[0][0] * gamma_1 + state[0][1] * gamma_2).tolist(), state[1]] for state in spectrum[1:-1]] ### Dropping the first and last state, because they correspond to the parents
    elif m < 0:
        spectrum = KS2(-m,omega_2,omega_1)
        return list(reversed(spectrum[1:-1]))
        ### the following command would return the new states in the global basis, as opposed to the parent's basis
        # return [[(state[0][0] * gamma_2 + state[0][1] * gamma_1).tolist(), state[1]] for state in spectrum[1:-1]] ### Dropping the first and last state, because they correspond to the parents


def matrix_multiply(a,b):
    zip_b = zip(*b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def t_expand(expr):
    """
    Takes an expression f(x,y) and computes the Taylor expansion in x and y 
    up to bi-degree fixed by max_deg.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    x + x * y
    up to bi-degree 2.
    """
    return series(expr.subs([(x, t * x), (y, t * y)]), t, 0, max_deg).removeO().subs(t, 1)


def t_degree(expr):
    """
    Takes an expression f(x,y) and computes the Taylor expansion in x and y 
    up to bi-degree fixed by max_deg. Then, picks the term of highest bi-degree,
    and returns the value of the degree.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    t * x + t**2 * x * y
    up to bi-degree 2. Will return 2.
    """
    return degree(series(expr.subs([(x, t * x), (y, t * y)]), t, 0, max_deg).removeO(), t)


def t_min_degree(expr):
    """
    Takes an expression f(x,y) and computes the Taylor expansion in x and y 
    up to bi-degree fixed by max_deg. Then, picks the term of lowest bi-degree.
    For example: 
    x / (1 - y) 
    with value of max_deg = 3 will give 
    t * x + t**2 * x * y
    up to bi-degree 2. Will return x.
    """
    return LT((t**1000 * series(expr.subs([(x, x / t), (y, y / t)]), t, 0, max_deg).removeO()).expand(), t).subs(t, 1)


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
    # b_set = set(map(tuple,a))  #need to convert the inner lists to tuples so they are hashable
    # b = map(list,b_set) #Now convert tuples back into lists (maybe unnecessary?)
    return b


def X(gamma):
    return (x ** gamma[0]) * (y ** gamma[1])


def K(m, gamma, omega, expr):
    """
    This implements the KS operator, in the active form.
    The charges are understood to be expressed in a basis where 
    gamma_1 = [1,0] and gamma_2 = [0,1]
    the actual pairing matrix will be specified by a single integer,
    and reads
    [[0, m], [-m, 0]] 
    """
    x_prime = x * ((1 + (-1)**(m*gamma[0]*gamma[1]-1) * X(gamma)) ** (m * gamma[1] * omega))
    y_prime = y * ((1 + (-1)**(m*gamma[0]*gamma[1]-1) * X(gamma)) ** (- m * gamma[0] * omega))
    return t_expand(expr.subs([(x, x_prime), (y, y_prime)], simultaneous=True))

def S(m, data, expr):
    """
    This implements a sequence of KS operations.
    The formatting of "data" should be as follows:
    [ [gamma_1 , omega_1]  ,  [gamma_2 , omega_2] ] 
    phase ordered from right to left.
    """
    temp = expr
    for i in range(len(data)):
        gamma = data[-(i+1)][0]
        omega = data[-(i+1)][1]
        temp = K(m, gamma, omega, temp)
    return temp


def position_sublist(lst, sublst):
    for i in range(len(lst)):
        if sublst == lst[i]:
            return i

def KS2(m, omega_1, omega_2):    
    """
    Solves the basic KS identity: 
    K_(0,1)^{omega_2} K_(1,0)^{omega_2} = (???)
    with 
    < (1,0) , (0,1) > = m
    We use the procedure suggested in eq. (5.2) of WWC (wild wall crossing paper)
    """

    ### Start by checking if we already computed this
    case = position_sublist(stored_keys, [m, omega_1, omega_2])

    if not case == None:
        return stored_results[case]

    else:
        ### The reference data: the LHS of the KSWCF
        data = [ [[0,1], omega_2],   [[1,0],omega_1]]
        var_x = S(m, data, x)
        var_y = S(m, data, y)

        ### Finding the RHS
        ansatz = [ [[1,0],omega_1], [[0,1], omega_2]]

        while True:
            ans_x = S(m, ansatz, x)
            ans_y = S(m, ansatz, y)

            expr_x = t_min_degree(((var_x - ans_x) / x).expand()) 
            # print "expr_x = %s" % expr_x
            expr_y = t_min_degree(((var_y - ans_y) / y).expand()) 
            # print "expr_y = %s" % expr_y

            if expr_x == 0 and expr_y ==0:
                break

            else:
                all_monomials_x = monomials_of_degree(t_degree(expr_x))
                # print all_monomials_x
                all_coefficients_x = [Poly(expr_x).coeff_monomial(mono) for mono in all_monomials_x]
                # print all_coefficients_x
                all_charges_x = charges_of_degree(t_degree(expr_x))
                # print all_charges_x
                all_degeneracies_x = [-all_coefficients_x[i] / (m * (-1)**(m * all_charges_x[i][0] * all_charges_x[i][1]) * all_charges_x[i][1]) for i in range(len(all_coefficients_x))]
                # print all_degeneracies_x

                all_monomials_y = monomials_of_degree(t_degree(expr_y))
                # print all_monomials_y
                all_coefficients_y = [Poly(expr_y).coeff_monomial(mono) for mono in all_monomials_y]
                # print all_coefficients_y
                all_charges_y = charges_of_degree(t_degree(expr_y))
                # print all_charges_y
                all_degeneracies_y = [-all_coefficients_y[i] / (- m * (-1)**(m * all_charges_y[i][0] * all_charges_y[i][1]) * all_charges_y[i][0]) for i in range(len(all_coefficients_y))]
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

                ### Now add the newly found BPS states to the ansatz, and reorder according to charge slope
                ansatz = ansatz + all_states
                # print ansatz
                decorated = [[(state[0][1]+0.01)/(state[0][0]+0.01), state] for i, state in enumerate(ansatz)] ## decorating a list of states by the charge slope, adding +0.01 to avoid dividing by zero
                ansatz = [state for slope, state in sorted(decorated)]
                # print ansatz

        stored_keys.append([m, omega_1, omega_2])
        stored_results.append(ansatz)

        return ansatz


###### SOME TESTS #####

# print KS2(2,1,1)
# print KS2(3,1,1)

# print stored_keys
# print stored_results
# print KS2(2,1,1)

# data = [[[1,0], 1], [[-1,2], 1]] 
# print progeny_2(data)
# print type(progeny_2(data))

# data = [[[-1,2], 1], [[1,0], 1]] 
# print progeny_2(data)

# data = [[[1,0], 1], [[0,1], 1]] 
# print progeny_2(data)

# data = [[[1,0], 1], [[0,3], 1]] 
# print progeny_2(data)

# print KS2(3,1,1)

