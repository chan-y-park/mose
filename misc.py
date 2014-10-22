import logging
import numpy
import datetime
from scipy.optimize import fsolve
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