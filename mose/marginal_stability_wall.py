import logging
import ast
from misc import deep_reverse
from numpy import array

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

    def __init__(self, all_intersections, fibration):
        ### !!! CONCEPTUALLY WRONG !!! 
        ### An MS wall doenst have a single pair of charges, 
        ### but rather several pairs, related by monodromies!
        ### Improve this piece of data.!
        self.gauge_charges = all_intersections[0].gauge_charges 
        self.flavor_charges = all_intersections[0].flavor_charges 

        self.degeneracies = all_intersections[0].degeneracies
        self.genealogy = all_intersections[0].genealogy
        self.fibration = fibration
        self.points = all_intersections
        ### the following enhances the self.points attribute, 
        ### possibly by adding branch-points
        self.enhance_ms_wall()
        self.delete_duplicate_points()     
        ### Now we reorder the list of self.points,
        ### according to the actual shape of the wall
        self.reorganize()
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
        a_random_point = self.points[0]
        kwall_0, kwall_1 = a_random_point.parents

        if kwall_0.__class__.__name__ == 'PrimaryKWall' \
                    and kwall_1.__class__.__name__ == 'PrimaryKWall':
            branch_point_0 = kwall_0.parents[0]
            branch_point_1 = kwall_1.parents[0]
            self.points.append(branch_point_0)
            self.points.append(branch_point_1)
        
    def reorganize(self):
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

            ### !!!
            ### NEED TO FIX THIS NUMERICAL CONSTANT SOMEHOW 
            ### IF TOO LARGE, IT WILL PRODUCE THOSE SECANTS 
            ### IN THE INTERPOLATION OF INTERSECTION POINTS
            ### TO MAKE A MS-WALL.
            ### IF TOO SMALL, WE WILL LOSE POINTS ON THE MS-WALL.
            ### !!!
            search_range = 0.75 * max_distance

            self.semi_arc_1 = []
            self.semi_arc_2 = []

            seen = [leftmost_point.locus]
            
            ### Build the first semi-arc
            current_point = leftmost_point
            arc_is_finished = False
            while not arc_is_finished:
                found_new_closest_point = False
                epsilon = search_range
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
                    seen.append(closest_point.locus)
                    self.semi_arc_1.append(closest_point)
                    current_point = closest_point

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
                    seen.append(closest_point.locus)
                    self.semi_arc_2.append(closest_point)
                    current_point = closest_point

            reorganized_points = self.semi_arc_1[::-1] + [leftmost_point] \
                                                            + self.semi_arc_2
            if len(reorganized_points) == len(self.points):
                self.points = reorganized_points
                pass
            else:
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


def getkey_real(int_point):
    return int_point.locus.real


def build_ms_walls(k_wall_networks):
    """
    This function creates MS walls, by sifting through all the intersections.
    These are clustered accoding to a certain choice of data.
    """
    all_intersections = []
    fibration = k_wall_networks[0].fibration
    for kwn in k_wall_networks:
        all_intersections += kwn.intersections
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
        if not data[i] in seen:
            walls.append([all_intersections[i]]) #start a new wall
            seen.append(data[i])
        else:
            walls[seen.index(data[i])].append(all_intersections[i])


    return [MarginalStabilityWall(x, fibration) for x in walls]


