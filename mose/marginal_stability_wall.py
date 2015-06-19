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
        self.charges = all_intersections[0].charges 
        ### warning: self.charges is given in the format {'[-1, 2]', '[1, 0]'}
        self.degeneracies = all_intersections[0].degeneracies
        self.genealogy = all_intersections[0].genealogy
        self.fibration = fibration
        self.points = all_intersections
        ### the following enhances the self.points attribute, 
        ### possibly by adding branch-points
        
        ### !!! NOTE -- temoporarily switched off the enhancement !!!
        ### It is causing trouble and should be rethought very carefully!
        ###
        # self.enhance_ms_wall()     
        ###

        ### Now we reorder the list of self.points,
        ### according to the actual shape of the wall
        self.reorganize()
        self.locus = [point.locus for point in self.points]
        MarginalStabilityWall.count += 1

    def enhance_ms_wall(self):
        ### now we add branch-points if they belong to the MS-walls,
        ### this is determined by wether one of the parents of the MS-wall
        ### is a primary trajectory. We get this info from the genealogy.
        ### The parent will be a primary kwall if its genealogy
        ### is directly a string, namely the identifier of the 
        ### branch point from which it sources
        gen = self.genealogy
        if type(gen[0]) == str:
            self.points.append(branch_point_with_identifier(\
                                                    self.fibration, gen[0]))
        if type(gen[1]) == str:
            self.points.append(branch_point_with_identifier(\
                                                    self.fibration, gen[1]))

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
            search_range = 1.1 * max_distance

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
                            # print " %s this is not the leftmost point" % pt.locus
                            if (distance <= abs(pt.locus \
                                                    - leftmost_point.locus)):
                                # print "but the distance to current_pt is less than to the leftmost_pt"
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
                raise ValueError('Lost some of the MS Wall points '\
                                    + 'while reorganizing it.')


def branch_point_with_identifier(fibration, label):
    found = False
    for b_pt in fibration.branch_points:
        if b_pt.genealogy == label:
            found = True
            return b_pt
    if found == False:
        raise ValueError('Problem when enhancing an MS wall: '\
                        +'cannot find the branch point with label {}'\
                        .format(label))



def getkey_real(int_point):
    return int_point.locus.real


def opposite(data):
    ### recall that charges of an intersection point are formatted as 
    ### {'[1,0]','[0,-2]'}
    ### to turn the above into a format such as
    ### [[1,0],[0,-2]]
    ### we do the following
    charges = map(ast.literal_eval, list(data[0]))
    genealogy = data[1]
    ### formatting back to the previous convention
    opp_charges =  set(map(str, map(list, map(lambda x: -x, \
                                                    map(array, charges)))))
    ### apparently, no need to take the "mirror" genealogy: makes sense!
    # opp_genealogy = set(map(str, deep_reverse(genealogy)))
    opp_genealogy = genealogy
    return [opp_charges, opp_genealogy]


def build_ms_walls(k_wall_networks):
    """
    This function creates MS walls, by sifting through all the intersections.
    These are clustered accoding to the genealogies and their charges.
    """
    all_intersections = []
    fibration = k_wall_networks[0].fibration
    for kwn in k_wall_networks:
        all_intersections += kwn.intersections
    ### to distinguish wall types, use certain data, defined in the following
    data = [[x.charges, x.genealogy] for x in all_intersections]
    seen = []
    walls = []

    logging.info(
                    '-----------------\n'
                    'Building MS walls\n'
                    '-----------------'
                )

    for i in range(len(data)):
        if not (data[i] in seen or opposite(data[i]) in seen):
            walls.append([all_intersections[i]]) #start a new wall
            seen.append(data[i])
        elif opposite(data[i]) in seen:
            walls[seen.index(opposite(data[i]))].append(all_intersections[i])
        else:
            walls[seen.index(data[i])].append(all_intersections[i])


    return [MarginalStabilityWall(x, fibration) for x in walls]


