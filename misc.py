import logging
import numpy
import datetime

def complexify(y):
    """ complexifies an array of two reals """
    return y[0] + 1j * y[1]

def dsz_pairing(gamma_1, gamma_2, dsz_matrix):
    return numpy.dot(numpy.dot(gamma_1, dsz_matrix), gamma_2)

def formatted_date_time():
	today = datetime.date.today()
	now = datetime.datetime.now().time().strftime("%H.%M")
	return str(today) + '-' + str(now)
