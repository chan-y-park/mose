"""
Find intersections of K-walls and package each of them into 
IntersectionPoint class.

Uses general-purpose module, intersection.py
"""
import platform
import os
import numpy
import logging
import ctypes
import pdb
from itertools import chain
from misc import sort_parent_kwalls, id_generator
from misc import dsz_pairing
from genealogy import build_genealogy_tree
from intersection import (NoIntersection, find_intersection_of_segments,
                          remove_duplicate_intersection,)
from zones import build_charge_orbit, determine_zone

INTERSECTION_ACCURACY = 1e-1

class IntersectionPoint:
    """The IntersectionPoint class.

    Attributes: 
    - locus (point on moduli space), 
    - index_1 (index of intersection point within parent number 1, 
      important to determine the charge at intersection); 
    - index_2 (index of intersection point within parent number 2, 
      important to determine the charge at intersection); 
    - genealogy

    Arguments: 
    data (as a triplet of [u, index_1, index_2]), 
    parents (as list of trajectories, ie objects)    
    """

    def __init__(self, data, parents):

        self.locus = data[0]
        index_1 = data[1]
        index_2 = data[2]

        ### Sort the parents and the corresponding indices
        ### according to how they come to the intersection
        ### as (kwall_l, kwall_r) ~ (x, y)
        self.parents, self.indices = sort_parent_kwalls(
            parents, [index_1, index_2]
        )

        self.fibration = self.parents[0].fibration

        ### This check has been deferred to the diagnostics module
        ### and will only be performed upon specific request.
        # check_marginal_stability_condition(self)

        self.degeneracies = [self.parents[0].degeneracy,
                             self.parents[1].degeneracy]
        self.genealogy = build_genealogy_tree(self)
        self.phase = self.parents[0].phase
        self.identifier = id_generator()


        self.gauge_charges = [self.parents[0].charge(self.indices[0]),
                                    self.parents[1].charge(self.indices[1])]
        self.flavor_charges = [self.parents[0].flavor_charge(self.indices[0]),
                                self.parents[1].flavor_charge(self.indices[1])]

        self.zone = determine_zone(self.locus, self.fibration.zones)
        self.charge_orbit = build_charge_orbit(self, self.fibration)



    def __eq__(self, other):
        if (self.parents == other.parents and
            self.indices[0] == other.indices[0] and
            self.indices[1] == other.indices[1]):
            return True
        else:
            return False


def get_k_wall_turning_points(k_wall):
    """
    Return a list of indices of turning points of a K-Wall,
    i.e. dx/dy = 0 or dy/dx = 0.
    """
    tps = []

    xs = k_wall.get_xs()
    ys = k_wall.get_ys()

    if len(xs) < 3:
        return tps

    for t in range(1, len(xs)-1):
        if ((xs[t] - xs[t-1]) * (xs[t+1] - xs[t]) < 0 or
            (ys[t] - ys[t-1]) * (ys[t+1] - ys[t]) < 0):
            tps.append(t)
        
    return tps

def find_new_intersection_points(
    prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix
):
    try:
        linux_distribution = platform.linux_distribution()[0]
        if linux_distribution != '':
            return find_new_intersection_points_using_cgal(
                prev_k_walls, new_k_walls, prev_intersection_points, 
                dsz_matrix, linux_distribution=linux_distribution,
            )
        else:
            raise OSError
    except OSError:
        logging.warning('CGAL not available; switch from '
                        'find_new_intersection_points_using_cgal() to '
                        'find_new_intersection_points_using_interpolation().')
        return find_new_intersection_points_using_interpolation(
            prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix
        )

#def find_new_intersection_points(
#    prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix
#):
#    return find_new_intersection_points_using_interpolation(
#        prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix
#    )


def find_new_intersection_points_using_cgal(
    prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix,
    linux_distribution=None,
):
    """
    Find new wall-wall intersections using CGAL 2d curve intersection.
    This function compares each K-walls pairwise, thereby having
    O(n^2) performance in theory.
    """

    new_intersection_points = []

    lib_name = 'libcgal_intersection'
    #if linux_distribution == 'Ubuntu':
    #    lib_name += '_ubuntu'
    #elif linux_distribution == 'Scientific Linux':
    #    lib_name += '_het-math2'
    #else:
    #    raise OSError

    logging.info('Using CGAL to find intersections.')

    # Load CGAL shared library.
    libcgal_intersection = numpy.ctypeslib.load_library(
        lib_name, 
        os.path.dirname(os.path.realpath(__file__)) + '/cgal_intersection/'
    )
    cgal_find_intersections_of_curves = (libcgal_intersection.
                                         find_intersections_of_curves)
    # Prepare types for CGAL library.
    array_2d_float = numpy.ctypeslib.ndpointer(
        dtype=numpy.float64,
        ndim=2,
        flags=['C_CONTIGUOUS', 'ALIGNED'],
    )

    cgal_find_intersections_of_curves.restype = ctypes.c_int
    cgal_find_intersections_of_curves.argtypes = [
        array_2d_float, ctypes.c_long,
        array_2d_float, ctypes.c_long,
        array_2d_float, ctypes.c_int,
    ]

    for n, k_wall_1 in enumerate(new_k_walls):
        # KWall.coordinates is a N-by-2 numpy array with dtype=float.
        curve_1 = numpy.require(k_wall_1.coordinates,
                                numpy.float64, ['A', 'C'])
        for k_wall_2 in chain(prev_k_walls, new_k_walls[n+1:]):
            if (
                # XXX: need to check local charges instead of initial ones.
                k_wall_1.parents == k_wall_2.parents or 
                k_wall_1 in k_wall_2.parents or
                k_wall_2 in k_wall_1.parents
            ):
                continue

            curve_2 = numpy.require(k_wall_2.coordinates,
                                    numpy.float64, ['A', 'C'])

            buffer_size = 10
            intersection_search_finished = False
            while not intersection_search_finished:
                intersections = numpy.empty((buffer_size, 2),
                                            dtype=numpy.float64)
                num_of_intersections = cgal_find_intersections_of_curves(
                    curve_1, ctypes.c_long(len(curve_1)),
                    curve_2, ctypes.c_long(len(curve_2)),
                    intersections, buffer_size,
                )
                if num_of_intersections > buffer_size:
                    logging.info('Number of intersections larger than '
                                 'the buffer size; increase its size '
                                 'and run intersection finding again.')
                    buffer_size = num_of_intersections
                else:
                    intersections.resize((num_of_intersections,2))
                    intersection_search_finished = True
            for ip_x, ip_y in intersections:
                # Find where to put the intersection point in the
                # given segment. It should be put AFTER the index
                # found below.
                curve_1_xs = k_wall_1.get_xs()
                x_full_1 = numpy.full(len(curve_1_xs), ip_x)
                i_a, i_b = numpy.argsort(abs(curve_1_xs - x_full_1))[:2]
                if i_a < i_b:
                    index_1 = i_a
                else:
                    index_1 = i_b
               
                curve_2_xs = k_wall_2.get_xs()
                x_full_2 = numpy.full(len(curve_2_xs), ip_x)
                i_c, i_d = numpy.argsort(abs(curve_2_xs - x_full_2))[:2]
                if i_c < i_d:
                    index_2 = i_c
                else:
                    index_2 = i_d

                logging.debug('intersection point: (%.8f, %.8f) '
                                'at index_1 = %d, index_2 = %d',
                                ip_x, ip_y, index_1, index_2)

                # Check local charges at the intersection.
                # XXX: What if the intersection happens 
                # at a location very close to a branch cut?
                if dsz_pairing(k_wall_1.charge(index_1),
                               k_wall_2.charge(index_2), 
                               dsz_matrix,) == 0: 
                    continue
                new_intersection_points.append(
                    IntersectionPoint(
                        [complex(ip_x + 1j*ip_y), index_1, index_2],
                        [k_wall_1, k_wall_2],
                    )
                )

    return new_intersection_points
    

def find_new_intersection_points_using_interpolation(
    prev_k_walls, new_k_walls, prev_intersection_points, dsz_matrix
):
    """
    Find new wall-wall intersections using interpolations.

    This first gets segments of a K-wall by cutting it at its
    turning points, and then find a SINGLE intersection between
    two segments of K-walls. Therefore it will fail to find more
    than one intersection per pair of segments.
    """
    logging.info('Using NumPy interpolation to find intersections.')

    new_intersection_points = []

    for n, k_wall_1 in enumerate(new_k_walls):
        tps_1 = get_k_wall_turning_points(k_wall_1)
        xs_1 = k_wall_1.get_xs()
        ys_1 = k_wall_1.get_ys()
        split_1 = [0] + tps_1 + [len(xs_1)-1]

        for k_wall_2 in chain(prev_k_walls, new_k_walls[n+1:]):
            # Exclude some cases from being checked for 
            # intersections, see the if statements below.
            if (
                # XXX: need to check local charges instead of initial ones.
                #dsz_pairing(k_wall_1.charge(0), k_wall_2.charge(0), 
                #            dsz_matrix) == 0 or 
                k_wall_1.parents == k_wall_2.parents or 
                k_wall_1 in k_wall_2.parents or
                k_wall_2 in k_wall_1.parents
            ):
                continue

            tps_2 = get_k_wall_turning_points(k_wall_2)
            xs_2 = k_wall_2.get_xs()
            ys_2 = k_wall_2.get_ys()
            split_2 = [0] + tps_2 + [len(xs_2)-1]

            for i_1 in range(1, len(split_1)):
                t_1_i = split_1[i_1-1]
                t_1_f = split_1[i_1]
                seg_1_xs = xs_1[t_1_i:t_1_f+1] 
                seg_1_ys = ys_1[t_1_i:t_1_f+1]
                seg_1 = (seg_1_xs, seg_1_ys)

                for i_2 in range(1, len(split_2)):
                    t_2_i = split_2[i_2-1]
                    t_2_f = split_2[i_2]
                    seg_2_xs = xs_2[t_2_i:t_2_f+1] 
                    seg_2_ys = ys_2[t_2_i:t_2_f+1]
                    seg_2 = (seg_2_xs, seg_2_ys)

                    try:
                        ip_x, ip_y = find_intersection_of_segments(
                            seg_1, seg_2, accuracy=INTERSECTION_ACCURACY,
                        )
                        # Find where to put the intersection point in the
                        # given segment. It should be put AFTER the index
                        # found below.
                        x_full_1 = numpy.full((t_1_f+1-t_1_i,), ip_x)
                        i_a, i_b = numpy.argsort(abs(seg_1_xs - x_full_1))[:2]
                        if i_a < i_b:
                            index_1 = t_1_i + i_a
                        else:
                            index_1 = t_1_i + i_b
                       
                        x_full_2 = numpy.full((t_2_f+1-t_2_i,), ip_x)
                        i_c, i_d = numpy.argsort(abs(seg_2_xs - x_full_2))[:2]
                        if i_c < i_d:
                            index_2 = t_2_i + i_c
                        else:
                            index_2 = t_2_i + i_d

                        logging.debug('intersection point: (%.8f, %.8f) '
                                        'at index_1 = %d, index_2 = %d',
                                        ip_x, ip_y, index_1, index_2)

                        # Check local charges at the intersection.
                        # XXX: What if the intersection happens 
                        # at a location very close to a branch cut?
                        if dsz_pairing(k_wall_1.charge(index_1),
                                       k_wall_2.charge(index_2), 
                                       dsz_matrix,) == 0: 
                            continue
                        new_intersection_points.append(
                            IntersectionPoint(
                                [complex(ip_x + 1j*ip_y), index_1, index_2],
                                [k_wall_1, k_wall_2],
                            )
                        )
                    except NoIntersection:
                        pass

    remove_duplicate_intersection(new_intersection_points, 
                                  prev_intersection_points)

    logging.info('Evaluating intersections of NEW Kwalls with ALL Kwalls: '
                 'found %d of them.', len(new_intersection_points))

    return new_intersection_points


# XXX: Deprecated
def find_new_intersections_using_hit_table(
    kwalls, new_kwalls, intersections, hit_table, dsz_matrix
):
    """
    Find new wall-wall intersections using HitTable.
    """
    logging.info('Using HitTable & NumPy interpolation to find intersections.')

    new_ints = []
    i_0 = len(kwalls)

    for i, traj in list(enumerate(new_kwalls)):
        # In hit_table, curves with an index higher than i_0 are
        # the new K-walls being added here.
        hit_table.fill(i_0+i, traj.coordinates)

    for bin_key in hit_table:

        if len(hit_table[bin_key]) == 1:
            # Only one curve in the bin. Skip the rest.
            continue
        # More than two curves hit the bin.
        logging.debug('bin_key with hit: %s', bin_key)
        for i_1, i_2 in combinations(hit_table[bin_key], 2):
            # logging.debug('i_0, i_1, i_2 = %s, %s, %s', i_0, i_1, i_2)
            # NOTE Chan: to get self-intersection, use 
            # combinations_with_replacement.
            if i_1 < i_0 and i_2 < i_0:
                # We checked this intersection already. Skip the rest.
                continue
            # At least one of the KWall is a new one.   
            if i_1 >= i_0:
                kwall_1 = new_kwalls[i_1 - i_0]
            else:
                kwall_1 = kwalls[i_1]

            if i_2 >= i_0:
                kwall_2 = new_kwalls[i_2 - i_0]
            else:
                kwall_2 = kwalls[i_2]

            # I am excluding some cases from being checked for 
            # intersections, see the if statements below
            if (dsz_pairing(kwall_1.charge(0), kwall_2.charge(0), 
                            dsz_matrix) == 0 or 
                # NOTE: We have to worry about losing an intersection 
                # of two K-walls from the same intersection
                # by imposing the following condition, which implies
                # that there will not be an intersection more than once
                # inside a bin.
                kwall_1.parents == kwall_2.parents or 
                kwall_1 in kwall_2.parents or
                kwall_2 in kwall_1.parents):
                continue

            list_of_intersection_points = []

            for t1_i, t1_f in hit_table[bin_key][i_1]:
                segment_1 = kwall_1.coordinates[t1_i:t1_f+1].T
                for t2_i, t2_f in hit_table[bin_key][i_2]:
                    segment_2 = kwall_2.coordinates[t2_i:t2_f+1].T
                    logging.debug('t1_i, t1_f = %d, %d', t1_i, t1_f) 
                    logging.debug('t2_i, t2_f = %d, %d', t2_i, t2_f) 
                    try:
                    # NOTE: we assume there is only one intersection
                    # between two segments, hoping that we divided
                    # the curves as many times as required for the
                    # assumption to hold.
                        intersection_point = find_intersection_of_segments(
                            segment_1, segment_2,
                            hit_table.get_bin_location(bin_key),
                            hit_table._bin_size
                        )

                        ipx, ipy = intersection_point
                        # Find where to put the intersection point in the
                        # given segment. It should be put AFTER the index
                        # found below.
                        
                        dt_1 = min(
                            range(len(segment_1)), 
                            key=lambda i: abs(segment_1[i][0]-ipx)
                        )
                        logging.debug('dt_1 = %d', dt_1)
                        if dt_1-1 >=0 and (
                            (segment_1[dt_1, 0] < ipx < segment_1[dt_1-1, 0]) 
                            or
                            (segment_1[dt_1, 0] > ipx > segment_1[dt_1-1, 0])
                            ):
                            dt_1 -= 1
                        index_1 = t1_i + dt_1 

                        dt_2 = min(
                            range(len(segment_2)), 
                            key=lambda i: abs(segment_2[i][0]-ipx)
                        )
                        logging.debug('dt_2 = %d', dt_2)
                        if dt_2-1 >=0 and (
                            (segment_2[dt_2, 0] < ipx < segment_2[dt_2-1, 0]) 
                            or
                            (segment_2[dt_2, 0] > ipx > segment_2[dt_2-1, 0])
                            ):
                            dt_2 -= 1
                        index_2 = t2_i + dt_2 
                        
                        logging.debug('intersection point: (%.8f, %.8f) '
                                        'at index_1 = %d, index_2 = %d',
                                        ipx, ipy, index_1, index_2)

                        list_of_intersection_points.append(
                            [complex(ipx + 1j*ipy), index_1, index_2]
                        )
                    except NoIntersection:
                        pass


            new_ints += [IntersectionPoint(intersection, [kwall_1, kwall_2]) 
                         for intersection in list_of_intersection_points] 

    remove_duplicate_intersection(new_ints, intersections)

    logging.info('Evaluating intersections of NEW Kwalls with ALL Kwalls: '
                    'found %d of them.', len(new_ints))

    return new_ints
