"""
Find intersections of K-walls and package each of them into 
IntersectionPoint class.

Uses general-purpose module, intersection.py
"""
import logging
from itertools import combinations
from misc import dsz_pairing
from config import DSZ_MATRIX
from intersection import NoIntersection, find_intersection_of_segments
from genealogy import build_genealogy_tree

class IntersectionPoint:
    """
    The IntersectionPoint class.

    Attributes: 
        locus (point on moduli space), 
        index_1 (index of intersection point within parent number 1, 
                    important to determine the charge at intersection), 
        index_2 (index of intersection point within parent number 2, 
                    important to determine the charge at intersection), 
    genealogy

    Arguments: 
        data (as a triplet of [u, index_1, index_2]), 
        parents (as list of trajectories, ie objects)    
    """

    def __init__(self,data,parents):

        self.parents = parents
        self.locus = data[0]
        self.index_1 = data[1]
        self.index_2 = data[2]
        self.charges = {
            str(parents[0].charge(self.index_1)), 
            str(parents[1].charge(self.index_2))
        }
        ### note the { } and conversion to strings, 
        ### since the charges are useful for classification purposes, mostly
        self.degeneracies = [parents[0].degeneracy, parents[1].degeneracy]
        self.genealogy = build_genealogy_tree(self)
        self.phase = parents[0].phase

    def __eq__(self, other):
        if (self.parents == other.parents and
            self.index_1 == other.index_1 and
            self.index_2 == other.index_2):
            return True
        else:
            return False

def remove_duplicate_intersection(new_ilist, old_ilist):
    """
    Remove any duplicate in new_ilist1, then remove intersections 
    of new_ilist that already exist in old_ilist
    """
    temp_ilist = new_ilist 

    for intersection1, intersection2 in combinations(temp_ilist, 2):
        if intersection1 == intersection2:
            new_ilist.remove(intersection2)
        else:
            continue

    temp_ilist = new_ilist 

    for new_intersection in temp_ilist:
        for intersection in old_ilist:
            if new_intersection == intersection:
                new_ilist.remove(new_intersection)

def find_new_intersections(kwalls, new_kwalls, intersections, hit_table):

    new_ints = []
    i_0 = len(kwalls)

    for i, traj in list(enumerate(new_kwalls)):
        hit_table.fill(i_0+i, traj.coordinates)

    for bin_key in hit_table:

        if len(hit_table[bin_key]) == 1:
            # Only one curve in the bin. Skip the rest.
            continue
        # More than two curves hit the bin.
        logging.debug('bin_key with hit: %s', bin_key)
        for i_1, i_2 in combinations(hit_table[bin_key], 2):
            #logging.debug('i_0, i_1, i_2 = %s, %s, %s', i_0, i_1, i_2)
            # NOTE Chan: to get self-intersection, use 
            # combinations_with_replacement.
            if i_1 < i_0 and i_2 < i_0:
                # We checked this intersection already. Skip the rest.
                continue
            # At least one of the KWall is a new one.   
            if i_1 >= i_0:
                traj_1 = new_kwalls[i_1 - i_0]
            else:
                traj_1 = kwalls[i_1]

            if i_2 >= i_0:
                traj_2 = new_kwalls[i_2 - i_0]
            else:
                traj_2 = kwalls[i_2]

            # NOTE Pietro: I am excluding some cases from being checked for 
            # intersections, see the if statements below
            if (dsz_pairing(traj_1.charge(0), traj_2.charge(0), 
                    DSZ_MATRIX) == 0 or                 
                traj_1.parents == traj_2.parents or 
                traj_1 in traj_2.parents):
                continue

            list_of_intersection_points = []

            for t1_i, t1_f in hit_table[bin_key][i_1]:
                segment_1 = traj_1.coordinates[t1_i:t1_f+1]
                for t2_i, t2_f in hit_table[bin_key][i_2]:
                    segment_2 = traj_2.coordinates[t2_i:t2_f+1]
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
                            key=lambda i: abs(segment_1[i, 0]-ipx)
                        )
                        logging.debug('dt_1 = %d', dt_1)
                        if (dt_1-1 >=0) and (
                            (segment_1[dt_1, 0] < ipx < segment_1[dt_1-1, 0]) 
                            or
                            (segment_1[dt_1, 0] > ipx > segment_1[dt_1-1, 0])
                        ):
                            dt_1 -= 1
                        index_1 = t1_i + dt_1 

                        dt_2 = min(
                            range(len(segment_2)), 
                            key=lambda i: abs(segment_2[i, 0]-ipx)
                        )
                        logging.debug('dt_2 = %d', dt_2)
                        if (dt_2-1 >=0) and (
                            (segment_2[dt_2, 0] < ipx < segment_2[dt_2-1, 0]) 
                            or
                            (segment_2[dt_2, 0] > ipx > segment_2[dt_2-1, 0])
                        ):
                            dt_2 -= 1
                        index_2 = t2_i + dt_2 
                        
                        logging.debug('intersecion point: (%.8f, %.8f) '
                                        'at index_1 = %d, index_2 = %d',
                                        ipx, ipy, index_1, index_2)
                        list_of_intersection_points.append(
                            [ipx + 1j*ipy, index_1, index_2]
                        )
                    except NoIntersection:
                        pass



            new_ints += [
                IntersectionPoint(intersection, [traj_1, traj_2]) 
                for intersection in list_of_intersection_points
            ] 

    remove_duplicate_intersection(new_ints, intersections)

    logging.info('Evaluating intersections of NEW Kwalls with ALL Kwalls: '
                    'found %d of them.', len(new_ints))

    return new_ints
    


