import logging 
import ast 
from misc import deep_reverse 
from numpy import array
from zones import orbit_is_contained, unite_orbits 
from cmath import phase, pi


### IF TOO LARGE, IT WILL PRODUCE THOSE SECANTS 
### IN THE INTERPOLATION OF INTERSECTION POINTS
### TO MAKE A MS-WALL.
### IF TOO SMALL, WE WILL LOSE POINTS ON THE MS-WALL.
SEARCH_RANGE_RATIO = 0.15 

from parameters import (ENHANCE_MS_WALLS, SORT_BY_GENEALOGY, MS_WALLS_SORTING,
                        CLEANUP_MS_WALLS,)

class MarginalStabilityWall:
    """
    Th MS wall class.
    Attributes:
    charges, degeneracy -- those of (any) two K-walls meeting theta_range
    points -- the list of actual intersectionpoint objects
    locus -- the list of points (UNSORTED FOR NOW)
    genealogy -- that of any intersection making up the wall

    Note: intersections are sorted following the arc traced by the ms-wall, 
    moreover (depending on the case) MS walls are enhanced with singularities 
    on which they end. 
    Such singularities are added at the end of the argument all_intersections, 
    they are given as instances of the branch-point class, and are NOT 
    converted to intersection-point objects.
    """

    count = 0

    def __init__(self, all_intersections, fibration, charge_orbit):
        ### !!! CONCEPTUALLY WRONG !!! 
        ### An MS wall doenst have a single pair of charges, 
        ### but rather several pairs, related by monodromies!
        ### Improve this piece of data.!
        self.gauge_charges = all_intersections[0].gauge_charges 
        self.flavor_charges = all_intersections[0].flavor_charges 

        ### The charge orbit of an MS wall, like that of an
        ### intersection point, is a list of elements
        ### [... , [[g1, g2], [f1, f2]] , ...]
        ### where g1, g2 are the gauge charges of the parent kwalls
        ### while f1, f2 are their flavor charges.
        ### The list contains elements of an orbit of these charges
        ### by the SL(2,Z) subgroup generated by the monodromies.
        ### By construction, the charge orbit of an MS wall 
        ### is larger or equal to the orbits of all it intersection
        ### points.
        self.charge_orbit = charge_orbit

        self.degeneracies = all_intersections[0].degeneracies
        self.genealogy = all_intersections[0].genealogy
        self.fibration = fibration
        self.points = all_intersections
        self.delete_duplicate_points()     

        self.deleted_points = []

        if MS_WALLS_SORTING == 'sweep':
            ### Now we reorder the list of self.points,
            ### according to the actual shape of the wall
            self.reorganize_sweep_sorting()
        elif MS_WALLS_SORTING == 'neighbor':
            ### Now we reorder the list of self.points,
            ### according to the actual shape of the wall
            self.reorganize_by_nearest_neighbor()
        elif MS_WALLS_SORTING == 'phase':
            ### Now we reorder the list of self.points,
            ### according to the phase \zeta
            self.reorganize_phase_sorting()

        if ENHANCE_MS_WALLS == True:
            ### the following enhances the self.points attribute, 
            ### possibly by adding branch-points
            self.enhance_ms_wall()
        else:
            pass
        
        if CLEANUP_MS_WALLS == True:
            self.cleanup()
        else:
            pass

        self.locus = [point.locus for point in self.points]
        MarginalStabilityWall.count += 1

    def delete_duplicate_points(self):
        seen = set()
        uniq = []
        for x in self.points:
            if x.locus not in seen:
                uniq.append(x)
                seen.add(x.locus)
        self.points = uniq

    def enhance_ms_wall(self):
        ### now we add branch-points if they belong to the MS-walls,
        ### this is determined by wether BOTH the parents of the MS-wall
        ### are primary kwalls.
        ### This function is called AFTER the points of the MS wall have 
        ### been sorted
        first_point = self.points[0]
        kwall_0, kwall_1 = first_point.parents

        if kwall_0.__class__.__name__ == 'PrimaryKWall' \
                    and kwall_1.__class__.__name__ == 'PrimaryKWall':
            branch_point_0 = kwall_0.parents[0]
            branch_point_1 = kwall_1.parents[0]
            
            dist_0_0 = abs(branch_point_0.locus - self.points[0].locus)
            dist_0_last = abs(branch_point_0.locus - self.points[-1].locus)
            if dist_0_0 < dist_0_last:
                self.points.insert(0, branch_point_0)
                self.points.append(branch_point_1)
            else:
                self.points.insert(0, branch_point_1)
                self.points.append(branch_point_0)

    
    ### OLD METHOD: This starts from the leftmost point (smallest Real part)
    ### of the intersection points of the MS wall, then tracks the 
    ### upper and lower arcs separately, by looking repeatedly for the 
    ### nearest neighbor, within a certain radius that we specify.
    def reorganize_by_nearest_neighbor(self):
        """
        MS walls are arcs, this function reorders the intersections, 
        following the arc.
        """
        if len(self.points) <= 2:
            pass

        else:
            points_copy = [pt for pt in self.points]
            leftmost_point = sorted(points_copy, key=getkey_real)[0]
            max_distance = max([abs(pt_1.locus - pt_2.locus) \
                                                for pt_1 in self.points \
                                                for pt_2 in self.points])    
            min_distance = min([abs(pt_1.locus - pt_2.locus) \
                                                for pt_1 in self.points \
                                                for pt_2 in self.points \
                                                if pt_1.locus != pt_2.locus])    

            # search_range = SEARCH_RANGE_RATIO * max_distance

            ### This is just a starting value.
            ### The search range will then be dynamically adapted.
            search_range = 1.1 * max_distance

            self.semi_arc_1 = []
            self.semi_arc_2 = []

            seen = [leftmost_point.locus]
            
            ### Build the first semi-arc
            current_point = leftmost_point
            arc_is_finished = False
            while not arc_is_finished:
                found_new_closest_point = False
                
                ### Set the maximal distance of points to be considered
                ### in a neighborhood of current_point
                epsilon = search_range
                
                ### We scan over all points for the next closest one
                for pt in self.points:
                    distance = abs(pt.locus - current_point.locus)
                    if distance < epsilon and (not (pt.locus in seen)): 
                        ### Now we check that when we actually
                        ### get at the end of the semi-arc 1, we don't 
                        ### switch to the beginning of the semi-arc 2.
                        ### The check involves comparing the distance
                        ### pt - current_point
                        ### with
                        ### current_point - leftmost_locus
                        ### The former should be smaller than the latter
                        ### However, we must skip this check at the 
                        ### very beginning, because in that case
                        ### current_point = leftmost_locus!
                        if current_point.locus != leftmost_point.locus:
                            if (distance <= abs(pt.locus \
                                                    - leftmost_point.locus)):
                                found_new_closest_point = True
                                epsilon = distance
                                closest_point = pt
                            else:
                                pass
                        else:
                            found_new_closest_point = True
                            epsilon = distance
                            closest_point = pt

                if found_new_closest_point == False:
                    arc_is_finished = True
                else:
                    ### We update the search range for the next iteration
                    search_range = 2.0 * abs(closest_point.locus - \
                                                        current_point.locus)
                    seen.append(closest_point.locus)
                    self.semi_arc_1.append(closest_point)
                    current_point = closest_point


            ### For the second semi-arc, the appropriate search range
            ### to begin with is (twice) the distance between the leftmost
            ### point and its direct neighbor in arc_1
            search_range = 2.0 * abs(leftmost_point.locus - \
                                                    self.semi_arc_1[0].locus)

            ### Build the second semi-arc
            current_point = leftmost_point
            arc_is_finished = False
            while not arc_is_finished:
                found_new_closest_point = False
                epsilon = search_range
                for pt in self.points:
                    distance = abs(pt.locus - current_point.locus)
                    if distance < epsilon and not (pt.locus in seen):
                        found_new_closest_point = True
                        epsilon = distance
                        closest_point = pt

                if found_new_closest_point == False:
                    arc_is_finished = True
                else:
                    ### We update the search range for the next iteration
                    search_range = 2.0 * abs(closest_point.locus - \
                                                        current_point.locus)
                    seen.append(closest_point.locus)
                    self.semi_arc_2.append(closest_point)
                    current_point = closest_point

            reorganized_points = self.semi_arc_1[::-1] + [leftmost_point] \
                                                            + self.semi_arc_2
            if len(reorganized_points) == len(self.points):
                self.points = reorganized_points
                pass
            else:
                self.points = reorganized_points
                print "\n IN WALL #%s" % self.count
                print "\noriginal points"
                # print self.points
                print [x.locus for x in self.points]
                print "\nreorganized points"
                # print reorganized_points
                print [x.locus for x in reorganized_points]
                logging.info('Lost some of the MS Wall points '\
                                    + 'while reorganizing it.\n' \
                                    + 'Try increasing the search radius.')
    

    ### NEW METHOD: an MS wall is assumed to have the shape of an arc, 
    ### we fix a point somewhere in the middle of it, and sort 
    ### all the intersection points accordint to the phase of the 
    ### radius connecting them to the basepoint.
    ### Then, we figure out where are the endpoints and reorganize 
    ### accordingly.

    def reorganize_sweep_sorting(self):
        """
        MS walls are arcs, this function reorders the intersections, 
        following the arc.
        """
        if len(self.points) <= 2:
            pass

        else:
            basepoint = sum([pt.locus for pt in self.points]) \
                                                            / len(self.points)
            
            sweep_sorted_pts = [x[0] for x in sorted(
                                    [[pt, phase(pt.locus - basepoint)] \
                                    for pt in self.points], key=getkey_second)
                                ]
            
            gap_start, gap_end = find_gap_for_sweep_sorting(sweep_sorted_pts, \
                                                                    basepoint)

            

            ### If the gap lies across the negative real axis, then 
            ### our points are already sorted correctly
            if gap_start == len(sweep_sorted_pts) and gap_end == 0:
                self.points = sweep_sorted_pts

            ### Otherwise we need to cut-and-paste
            else:
                ### recall that the phase is defined between -pi and pi,
                ### so we should cut and paste accordingly
                reorganized_points = sweep_sorted_pts[gap_end : ] \
                                            + sweep_sorted_pts[ : gap_start+1]

                if len(reorganized_points) == len(self.points):
                    self.points = reorganized_points
                    pass
                else:
                    logging.info("\n In MS wall #{} lost some of the points '\
                                + 'while reorganizing it.'".format(self.count))
                    logging.info("\noriginal points\n{}"\
                                    .format([x.locus for x in self.points]))
                    logging.info("\nreorganized points"\
                                .format([x.locus for x in reorganized_points])) 
                    self.points = reorganized_points

    ### NEW METHOD: an MS wall is assumed to have the shape of an arc, 
    ### we fix a point somewhere in the middle of it, and sort 
    ### all the intersection points accordint to the phase of the 
    ### radius connecting them to the basepoint.
    ### Then, we figure out where are the endpoints and reorganize 
    ### accordingly.

    def reorganize_phase_sorting(self):
        """
        Intersections come with a phase, that of the kwall network they belong 
        to. This algorithm uses that phase to order intersections by phase.
        """
        if len(self.points) <= 2:
            pass

        else:
            points_copy = [pt for pt in self.points]
            phase_sorted_pts = sorted(points_copy, key=getkey_phase)
            
            gap_start, gap_end = find_gap_for_phase_sorting(phase_sorted_pts)

            

            ### If the gap lies across the negative real axis, then 
            ### our points are already sorted correctly
            if gap_start == len(phase_sorted_pts) and gap_end == 0:
                self.points = phase_sorted_pts

            ### Otherwise we need to cut-and-paste
            else:
                ### recall that the phase is defined between -pi and pi,
                ### so we should cut and paste accordingly
                reorganized_points = phase_sorted_pts[gap_end : ] \
                                            + phase_sorted_pts[ : gap_start+1]

                if len(reorganized_points) == len(self.points):
                    self.points = reorganized_points
                    pass
                else:
                    logging.info("\n In MS wall #{} lost some of the points '\
                                + 'while reorganizing it.'".format(self.count))
                    logging.info("\noriginal points\n{}"\
                                    .format([x.locus for x in self.points]))
                    logging.info("\nreorganized points"\
                                .format([x.locus for x in reorganized_points])) 
                    self.points = reorganized_points

    def cleanup(self):
        ### Must be called after the wall is sorted. Preferably
        ### sweep-sorted, or phase-sorted.
        ### For each point, figures if its distance from its 
        ### neighbors is larger than the distance bewteen neighbors
        ### themselves, in that case it drops the point.
        ### The points are labeled 1,2,3 for previous, currrent, following.
        
        ### SHOULD IMPROVE with a criterion for selecting the 1st and last point
        ### right now they are picked up by default.
        
        if len(self.points) < 5:
            pass

        else:
            select_points = []

            ### compute the average distance between consecutive points,
            ### excluding the endpoints (which could be branch-points)
            tot_avg_dist = sum([ \
                        abs(self.points[j].locus - self.points[j-1].locus) \
                        for j in range(2, len(self.points) - 1, 1)]) / \
                        (len(self.points) - 2)

            for i in range(1, len(self.points) - 1, 1):
                dist_12 = abs(self.points[i].locus - self.points[i-1].locus)
                dist_23 = abs(self.points[i].locus - self.points[i+1].locus)
                dist_13 = abs(self.points[i-1].locus - self.points[i+1].locus)

                if dist_13 < dist_12 or dist_13 < dist_23:
                    self.deleted_points.append(self.points[i])
                    pass
                elif dist_12 > 2.0 * tot_avg_dist or \
                                            dist_23 > 2.0 * tot_avg_dist:
                    ### if the distance from one of the neighboring points 
                    ### is larger than 2 times the average distance so far 
                    self.deleted_points.append(self.points[i])
                    pass
                else:
                    select_points.append(self.points[i])

            select_points.insert(0, self.points[0])
            select_points.append(self.points[-1])

            self.points = select_points

        pass


def find_gap_for_sweep_sorting(sweep_sorted_pts, basepoint):
    ### we compute first the phase gap between the
    ### points across the negative real axis
    pt_0 = sweep_sorted_pts[0]
    pt_1 = sweep_sorted_pts[-1]
    phase_gap = 2.0 * pi + phase(pt_0.locus - basepoint) \
                                                - phase(pt_1.locus - basepoint)
    gap_start = len(sweep_sorted_pts)
    gap_end = 0

    ### Then we study angle differences between consecutive 
    ### intersection points, and check if any of these
    ### is larger than the above phase-gap.
    angle = phase(pt_0.locus - basepoint)
    for i, pt in enumerate(sweep_sorted_pts):
        new_angle = phase(pt.locus - basepoint)
        delta_angle = new_angle - angle
        angle = new_angle
        if delta_angle > phase_gap:
            phase_gap = delta_angle
            gap_start = i - 1
            gap_end = i
    
    return gap_start, gap_end


def find_gap_for_phase_sorting(phase_sorted_pts):    
    min_phase = phase_sorted_pts[0].phase
    max_phase = phase_sorted_pts[-1].phase
    if max_phase - min_phase > pi:
        raise ValueError('Cannot handle organizing MS walls generated '+ \
                'by kwall networks of phase-range larger than pi')

    ### We first compute the phase gap between the
    ### first and the last point. Here we use the crucial assumption 
    ### that the phase-range is at most pi.
    pt_0 = phase_sorted_pts[0]
    pt_1 = phase_sorted_pts[-1]
    phase_gap = pi + pt_0.phase - pt_1.phase
    gap_start = len(phase_sorted_pts)
    gap_end = 0

    ### Then we study angle differences between consecutive 
    ### intersection points, and check if any of these
    ### is larger than the above phase-gap.
    angle = pt_0.phase
    for i, pt in enumerate(phase_sorted_pts):
        new_angle = pt.phase
        delta_angle = new_angle - angle
        angle = new_angle
        if delta_angle > phase_gap:
            phase_gap = delta_angle
            gap_start = i - 1
            gap_end = i
    
    return gap_start, gap_end


def getkey_real(int_point):
    return int_point.locus.real


def getkey_second(x):
    return x[1]

def getkey_phase(x):
    return x.phase


def build_ms_walls(k_wall_networks):
    """
    This function creates MS walls, by sifting through all the intersections.
    These are clustered accoding to a certain choice of data.
    """
    all_intersections = []
    fibration = k_wall_networks[0].fibration
    for kwn in k_wall_networks:
        all_intersections += kwn.intersections
    
    if SORT_BY_GENEALOGY == True:
        logging.info('Constructing MS walls using genealogy')
    ### OLD METHOD, USED THE GENEALOGY
        ### to distinguish wall types, use the genealogy data
        data = [x.genealogy for x in all_intersections]

        seen = []
        walls = []

        logging.info(
                        '-----------------\n'
                        'Building MS walls\n'
                        '-----------------'
                    )
        for i in range(len(data)):
            check = data[i] in seen
            if check == False:
                walls.append([all_intersections[i]]) #start a new wall
                seen.append(data[i])
            else:
                ### we add this intersection point to the corresponding
                ### marginal stability wall
                walls[seen.index(data[i])].append(all_intersections[i])
    

    else:
        ### NEW METHOD: CHARGE ORBITS
        data = [x.charge_orbit for x in all_intersections]

        seen = []
        walls = []

        logging.info(
                        '-----------------\n'
                        'Building MS walls\n'
                        '-----------------'
                    )
        for i in range(len(data)):
            i_th_charge_orbit = data[i]
            check = orbit_is_contained(seen, i_th_charge_orbit)
            if check == 'not contained':
                walls.append([all_intersections[i]]) #start a new wall
                seen.append(i_th_charge_orbit)
            else:
                index = check
                ### we possibly enlarge the reference 'seen' orbit 
                ### by taking the union with this one
                seen[index] = unite_orbits(seen[index], i_th_charge_orbit)
                ### then we add this intersection point to the corresponding
                ### marginal stability wall
                walls[index].append(all_intersections[i])

    ### MS walls need to have a charge orbit assigned to them
    ### the i-th 'seen' orbit corresponds by construction to the 
    ### i-th MS wall, and it is the union of the orbits
    ### of all its intersection points.
    return [MarginalStabilityWall(x, fibration, seen[i]) \
                                                for i, x in enumerate(walls)]


