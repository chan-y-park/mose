import logging

def complexify(y):
    """ complexifies an array of two reals """
    return y[0] + 1j * y[1]

def dsz_pairing(gamma_1, gamma_2, dsz_matrix):
    import numpy as np
    return np.dot(np.dot(gamma_1,dsz_matrix),gamma_2)


