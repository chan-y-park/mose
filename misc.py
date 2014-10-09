import logging
import numpy
import datetime
from sympy.abc import u


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