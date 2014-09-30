import logging

class MarginalStabilityWall:
    """
    Th MS wall class.
    Attributes:
    charges, degeneracy -- those of (any) two K-walls meeting theta_range
    points -- the list of actual intersectionpoint objects
    locus -- the list of points (UNSORTED FOR NOW)
    genealogy -- that of any intersection making up the wall

    Note: intersections are not sorted, moreover (depending on the case)
    MS walls are enhanced with singularities on which they end. 
    Such singularities are added at the end of the argument all_intersections, 
    they are given as instances of the branch-point class, and are NOT 
    converted to intersection-point objects.
    """

    count = 0

    def __init__(self, all_intersections):
        self.charges = all_intersections[0].charges 
        ### warning: self.charges is given in the format {'[-1, 2]', '[1, 0]'}
        self.degeneracies = all_intersections[0].degeneracies
        self.genealogy = all_intersections[0].genealogy
        self.points = all_intersections
        # the following enhances the self.points attribute, 
        # possibly by adding branch-points
        self.enhance_ms_wall()     
        self.locus = [point.locus for point in self.points]
        MarginalStabilityWall.count += 1

    def enhance_ms_wall(self):
        ### now we add branch-points if they belong to the MS-walls,
        ### this is determined by wether one of the parents of the MS-wall
        ### is a primary trajectory. We get this info from the genealogy
        gen = self.genealogy
        if gen[0].__class__.__name__ == 'BranchPoint':
            self.points.append(gen[0])
        if gen[1].__class__.__name__ == 'BranchPoint':
            self.points.append(gen[1])


def build_ms_walls(k_wall_networks):
    """
    This function creates MS walls, by sifting through all the intersections.
    These are clustered accoding to the genealogies and their charges.
    """
    all_intersections = []
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
        if not (data[i] in seen):
            walls.append([all_intersections[i]]) #start a new wall
            seen.append(data[i])
        else:
            walls[seen.index(data[i])].append(all_intersections[i])

    return [MarginalStabilityWall(x) for x in walls]


