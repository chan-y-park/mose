import cmath
from misc import complexify


def build_genealogy_tree(intersection):
    """
    this function will return the genealogy tree of an intersection \
    in the form of a list such as [[BP1,BP2],BP3] for an intersection with the\
    obvious history
    """
    
    parents = intersection.parents
    indices = intersection.indices

    ### determine who's mom and who's dad by relative orientation
    delta_z_1 = complexify(parents[0].coordinates[indices[0]+1]) - \
                                complexify(parents[0].coordinates[indices[0]])
    delta_z_2 = complexify(parents[1].coordinates[indices[1]+1]) - \
                                complexify(parents[1].coordinates[indices[1]])

    if cmath.phase(delta_z_1 / delta_z_2) > 0:
        dad = parents[0].initial_point
        mom = parents[1].initial_point
    else:
        dad = parents[1].initial_point
        mom = parents[0].initial_point

    return [dad.genealogy, mom.genealogy] 
