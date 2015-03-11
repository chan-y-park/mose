import cmath
from misc import complexify


def build_genealogy_tree(intersection):
    """
    this function will return the genealogy tree of an intersection \
    in the form of a list such as [[BP1,BP2],BP3] for an intersection with the\
    obvious history
    """
    
    parents = intersection.parents
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    
    ### determine who's mom and who's dad by relative orientation
    delta_z_1 = complexify(parents[0].coordinates[index_1+1]) - \
                                    complexify(parents[0].coordinates[index_1])
    delta_z_2 = complexify(parents[1].coordinates[index_2+1]) - \
                                    complexify(parents[1].coordinates[index_2])

    if cmath.phase(delta_z_1 / delta_z_2) > 0:
        dad = parents[0].initial_point
        mom = parents[1].initial_point
    else:
        dad = parents[1].initial_point
        mom = parents[0].initial_point

    return [dad.genealogy, mom.genealogy] 