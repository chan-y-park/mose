"""
Objects and functions to find intersections of real 1-dim curves 
on a real 2-dim plane.

This module is based on the following libraries:
NumPy 1.9.0
SciPy 0.14.0
SymPy 0.7.5
"""
from numpy import __version__ as numpy_version
from scipy import __version__ as scipy_version
from sympy import __version__ as sympy_version
from scipy.interpolate import interp1d
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq, newton
# NOTE:this module requires SymPy 0.7.5.
from sympy import Interval, Intersection

from math import ceil
from itertools import combinations
from warnings import warn

from structure import Trajectory, IntersectionPoint

# Library version checks.
if numpy_version < '1.9.0':
    message = (str(numpy_version) + ' is lower than 1.9.0; '  
                'this module may not work properly.')
    warn(message, Warning)

if scipy_version < '0.14.0':
    message = (str(scipy_version) + ' is lower than 0.14.0; '  
                'this module may not work properly.')
    warn(message, Warning)

if sympy_version < '0.7.5':
    message = (str(sympy_version) + ' is lower than 0.7.5; '  
                'this module may not work properly.')
    warn(message, Warning)

class NoIntersection(Exception):
    """
    An exception class to be raised when failed 
    to find an intersection.
    """
    def __init__(self, value = ''):
        self.value = value

    def __str__(self):
        return repr(self.value)

class HitTable:
    """
    Construct a hash table of a coarse-grained plane.

    The table is implemented as a Python dictionary, whose item 
    has a key of the form (x_ind, y_ind), i.e. (1, 2), and the 
    corresponding value is a bin of the plane, which in turn 
    is a dictionary of 
       key: value =  ind: [[t1_i, t1_f], [t2_i, t2_f], ...], 
    where c = list_of_curves[ind] gives a curve that hits the 
    bin, c[tn_i] is the starting coordinate of the nth segment 
    of the binned curve, and similarly c[tn_f] gives its endpoint.

    Parameters
    ----------
    plane_range : [[x_min, x_max], [y_min, y_max]]
    bin_size : size of one square bin
    bin_fill_offset: When set to 0, each segment goes over the 
        boundaries of bins, whereas a nonzero value corresponds to 
        trimming each end of segment by that amount. A nonzero offset 
        can result in missing an intersection when it is near the 
        boundary of a bin.

    Methods
    -------
    get_bin_key(location): returns the key of the bin at the location.
    get_bin_location(bin_key): returns the center location of the bin 
        corresponding to bin_key.
    put(bin_key, curve_index, segment_range): put a segment of a curve
        into the bin corresponding to bin_key.
    fill(list_of_curves): fill HitTable with the curves 
        in list_of_curves.
    """

    def __init__(self, plane_range=[[-1.0, 1.0], [-1.0, 1.0]], 
                    bin_size=1.0, bin_fill_offset=0):
        self._plane_range = plane_range
        self._bin_size = bin_size
        self._bin_fill_offset = bin_fill_offset
        self._hit_table = {}
        self._num_x_bin, self._num_y_bin = \
            map(lambda (r_min, r_max): int(ceil((r_max - r_min)/bin_size)), 
                plane_range)

    def __getitem__(self, key):
        return self._hit_table[key]

    def __iter__(self):
        return self._hit_table.iterkeys()

    def get_bin_key(self, location):
        x_coord, y_coord = location
        [[x_min, x_max], [y_min, y_max]] = self._plane_range

        bin_key_x, bin_key_y = \
            map(lambda coord, r_min: int((coord - r_min)/self._bin_size),
                location, (x_min, y_min))

        return (bin_key_x, bin_key_y)

    def get_bin_location(self, bin_key):
        [[x_min, x_max], [y_min, y_max]] = self._plane_range

        location = map(lambda coord, r_min: coord * self._bin_size + r_min + \
                        0.5*self._bin_size, bin_key, [x_min, y_min])
        
        return location

    def put(self, bin_key, curve_index, segment_range):

        if bin_key not in self._hit_table:
            bin_key_x, bin_key_y = bin_key
            # Check if bin_key is not a tuple of two integers.
            if not isinstance(bin_key_x, int) or not isinstance(bin_key_y, int):
                raise KeyError('bin_key ' + str(bin_key) + 
                                    ' not a tuple of two integers.')
            # Check if bin_key is within the search range.
            if not (0 <= bin_key_x < self._num_x_bin):
                raise KeyError('bin_key_x = ' + str(bin_key_x) + 
                                ' out of range: [0, ' + str(self._num_x_bin) + 
                                ']')
            if not (0 <= bin_key_y < self._num_y_bin):
                raise KeyError('bin_key_y = ' + str(bin_key_y) + 
                                ' out of range: [0, ' + str(self._num_y_bin) + 
                                ']')
            # bin_key is legitimate; create an empty dict of curve segments.
            self._hit_table[bin_key] = {}
        if curve_index not in self._hit_table[bin_key]:
            # Check if curve_index is an integer
            if not isinstance(curve_index, int):
                raise KeyError('curve_index ' + str(curve_index) + ' not an '
                                'integer.')
            # Create an empty list for segments
            self._hit_table[bin_key][curve_index] = []

        self._hit_table[bin_key][curve_index].append(segment_range)

    def fill(self, list_of_curves):
        for curve_index, curve in enumerate(list_of_curves):
            t_i = t_0 = 0
            x_0, y_0 = curve[t_0]
            prev_bin_key = self.get_bin_key([x_0, y_0])
            for t_n, [x_n, y_n] in enumerate(curve):
                bin_key = self.get_bin_key([x_n, y_n])
                # Cut the curve into a segment where it goes over the current 
                # bin or where it has a turning point.
                if prev_bin_key != bin_key:
                    # curve[t_n - 1] and curve[t_n] are in different bins. 
                    # NOTE: bin_fill_offset = 0 gives segments that go over the
                    # boundaries of bins. A nonzero offset can result in
                    # missing an intersection when it is near the boundary
                    # of a bin.
                    try:
                        t_f = t_n - self._bin_fill_offset
                        if t_i < t_f:
                            self.put(prev_bin_key, curve_index, [t_i, t_f])
                    except KeyError:
                        # [x_n, y_n] is outside the intersection search range.
                        print e.value
                        pass
                    t_i = t_n
                    prev_bin_key = bin_key
                elif is_turning_point(curve, t_n):
                    try:
                        self.put(prev_bin_key, curve_index, [t_i, t_n])
                    except KeyError as e:
                        # [x_n, y_n] is outside the intersection search range.
                        print e.value
                        pass
                    t_i = t_n
                else:
                    # curve[t_n] is not the end of the segment.
                    continue
            # Take care of the last segment, if any.
            if t_i < t_n:
                try:
                    self.put(prev_bin_key, curve_index, [t_i, t_n])
                except KeyError as e:
                    # [x_n, y_n] is outside the intersection search range.
                    print e.value
                    pass

def is_turning_point(curve, t):
    """Check whether curve[t] is the turning point or not."""
    t_max = len(curve) - 1
    if t_max < 2 or t == 0 or t == t_max:
        return False

    x, y = [list(c) for c in zip(*curve[t-1:t+2])] 

    # Check if dx/dy = 0
    if (x[1] - x[0]) * (x[2] - x[1]) < 0:
        return True
    # Check if dy/dx = 0
    elif (y[1] - y[0]) * (y[2] - y[1]) < 0:
        return True
    else:
        return False

def find_intersection_of_segments(segment_1, segment_2, bin_center, bin_size,
    newton_maxiter=5):
    """
    Find an intersection of two segments of curves in the same bin.
    
    First find interpolations of segments using scipy.interp1d and 
    use SciPy's Brent method to find an intersection. When this 
    doesn't work, use SciPy's polynomial interpolation and then use 
    secant method to find an intersection.

    Parameters
    ----------
    segment_1, segment_2: Segments to find their intersection. Each 
        segment is a list of [x, y].
    bin_center, bin_size: center location and the size of the bin 
        that contains the segments.
    newton_maxiter: Maximum number of iterations for secant method. 
        When increased, this gives a better accuracy of the 
        intersection but it also greatly reduces the performance 
        due to many cases when there is no intersection but 
        the module keeps trying to find one.
    """
    # First check if the two segments share any x- and y-range.
    x1_i, y1_i = segment_1[0]
    x1_f, y1_f = segment_1[-1]
    
    x2_i, y2_i = segment_2[0]
    x2_f, y2_f = segment_2[-1]

    bin_center_x, bin_center_y = bin_center

    x1_interval = Interval(*sorted([x1_i, x1_f]))
    x2_interval = Interval(*sorted([x2_i, x2_f]))
    bin_x_interval = Interval(bin_center_x - 0.5*bin_size,
                                bin_center_x + 0.5*bin_size)

    y1_interval = Interval(*sorted([y1_i, y1_f]))
    y2_interval = Interval(*sorted([y2_i, y2_f]))
    bin_y_interval = Interval(bin_center_y - 0.5*bin_size,
                                bin_center_y + 0.5*bin_size)

    x_range = Intersection(x1_interval, x2_interval, bin_x_interval)
    y_range = Intersection(y1_interval, y2_interval, bin_y_interval)

    if (x_range.is_EmptySet or y_range.is_EmptySet or x_range.is_FiniteSet or
        y_range.is_FiniteSet):
        # The segments and the bin do not share a domain and therefore
        # there is no intersection.
        raise NoIntersection() 

    f1 = interp1d(*zip(*segment_1))
    f2 = interp1d(*zip(*segment_2))
    delta_f12 = lambda x: f1(x) - f2(x)

    try:
        intersection_x = brentq(delta_f12, x_range.start, x_range.end)
        intersection_y = f1(intersection_x)
    except ValueError:
        """
        (f1 - f2) has the same sign at x_range.start & x_range.end
        use Newton's method instead, and for that purpose interpolate
        curves using polynomial interpolation.
        """
        # TODO: maybe spline interpolation is better because it can
        # provide derivatives of interpolatioin functions, but couldn't make 
        # it work yet.
        """
        # NOTE: cubic spline interpolation requires more than 3 points.
        if len(segment_1) <= 3 or len(segment_2) <= 3:
            # not enough data; stop finding an intersection
            raise NoIntersection
        # UnivariateSpline requires x to be an increasing array.
        segment_1.sort(0)
        segment_2.sort(0)
        f1 = UnivariateSpline(*zip(*segment_1))
        f1_prime = f1.derivative()
        f2 = UnivariateSpline(*zip(*segment_2))
        f2_prime = f2.derivative()
        x0 = 0.5*(x_range.start + x_range.end)

        delta_f12 = lambda x: f1(x) - f2(x)
        delta_f12_prime = lambda x: f1_prime(x) - f2_prime(x)

        intersection_x = newton(delta_f12, x0, delta_f12_prime(x0))
        """
        f1 = BarycentricInterpolator(*zip(*segment_1))
        f2 = BarycentricInterpolator(*zip(*segment_2))
        delta_f12 = lambda x: f1(x) - f2(x)

        x0 = 0.5*(x_range.start + x_range.end)

        try:
            intersection_x = newton(delta_f12, x0, maxiter = newton_maxiter)
        except RuntimeError:
            # Newton's method fails to converge; declare no intersection
            raise NoIntersection()

        # Check if the intersection is within the curve range. 
        # If not, the intersecion is not valid. 
        if intersection_x not in x_range:
            raise NoIntersection()
        intersection_y = f1(intersection_x)
        if intersection_y not in y_range:
            raise NoIntersection()

    return [intersection_x, float(intersection_y)]

def find_intersection(list_of_curves, hit_table, search_range, bin_size):
    """
    Find intersection points of two curves on a real 2-dim plane.
    Returns a list of intersection points
    """
    list_of_intersection_points = []

    for bin_key in hit_table:
        if len(hit_table[bin_key]) == 1:
            continue
        # more than two curves hit the bin.
        for curve_index_1, curve_index_2 in combinations(hit_table[bin_key], 2):
            #DEBUG: print 'bin @ ' + str(hit_table.get_bin_location(bin_key))
            # NOTE: to get self-intersection; use combinations_with_replacement.
            for t1_i, t1_f in hit_table[bin_key][curve_index_1]:
                segment_1 = list_of_curves[curve_index_1][t1_i:t1_f + 1]
                for t2_i, t2_f in hit_table[bin_key][curve_index_2]:
                    segment_2 = list_of_curves[curve_index_2][t2_i:t2_f + 1]
                    try:
                    # NOTE: we assume there is only one intersection
                    # between two segments, hoping that we divided
                    # the curves as many times as required for the
                    # assumption to hold.
                        list_of_intersection_points.append(
                            find_intersection_of_segments(
                                segment_1, segment_2,
                                hit_table.get_bin_location(bin_key),
                                bin_size
                            )
                        )
                    except NoIntersection:
                        pass

    return list_of_intersection_points
