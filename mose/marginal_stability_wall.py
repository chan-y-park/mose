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
        self.enhance_ms_wall()     
        ### Now we reorder the list of self.points,
        ### according to the actual shape of the wall
        self.reorganize()
        self.locus = [point.locus for point in self.points]
        MarginalStabilityWall.count += 1

    def enhance_ms_wall(self):
        ### now we add branch-points if they belong to the MS-walls,
        ### this is determined by wether one of the parents of the MS-wall
        ### is a primary trajectory. We get this info from the genealogy

        ### The genealogy is no longer built in this way!
        # gen = self.genealogy
        # if gen[0].__class__.__name__ == 'BranchPoint':
        #     self.points.append(gen[0])
        # if gen[1].__class__.__name__ == 'BranchPoint':
        #     self.points.append(gen[1])

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
        leftmost_point = sorted(self.points, key=getkey_real)[0]
        min_distance = min([abs(pt_1.locus - pt_2.locus) \
                                            for pt_1 in self.points \
                                            for pt_2 in self.points \
                                            if pt_1 != pt_2])
        max_distance = max([abs(pt_1.locus - pt_2.locus) \
                                            for pt_1 in self.points \
                                            for pt_2 in self.points \
                                            if pt_1 != pt_2])    
        search_range = 1.1 * max_distance
        semi_arc_1 = []
        semi_arc_2 = []

        seen = [leftmost_point.locus]
        
        ### Build the first semi-arc
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
                semi_arc_1.append(closest_point)

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
                semi_arc_2.append(closest_point)

        self.points = semi_arc_1[::-1] + [leftmost_point] + semi_arc_2
        pass


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


