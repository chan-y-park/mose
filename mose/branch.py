class BranchPoint:
    """The BranchPoint class.

    Attributes: locus, charge
    Arguments: locus, charge
    """

    count = 0

    def __init__(
        self, 
        locus=None,
        positive_charge=None,
        monodromy_matrix=None,
        positive_period=None,
        identifier=None,
    ):
        self.charge = positive_charge
        self.positive_period = positive_period
        self.locus = locus
        self.count = BranchPoint.count
        ### Old version: this would use the object itself 
        ### as the ancestor in the genealogy tree.
        ### It's a nice idea, but has some tension with
        ### multiprocessing.
        # self.genealogy = self
        ### New version: use a string as identifier.
        ### This should be more robust against issues
        ### arising from multiprocessing.
        self.genealogy = identifier
        
        self.monodromy_matrix = monodromy_matrix
        BranchPoint.count += 1

    def __str__(self):
        return 'Branch point info: charge %s, locus %s ' % \
            (self.charge, self.locus)


# class BranchCut:
#     """
#     The BranchCut class.
#     Attributes: locus, charge, monodromy_matrix
#     """
#     count = 0

#     def __init__(self, branch_point):
#         self.charge = branch_point.charge
#         self.monodromy_matrix = branch_point.monodromy_matrix
#         self.locus = branch_point.locus
#         BranchCut.count += 1


def minimum_distance(branch_points):
    loci = [bpt.locus for bpt in branch_points]
    min_dist = abs(loci[0]-loci[1])
    for i, z_i in enumerate(loci):
        for j, z_j in list(enumerate(loci))[i+1:]:
            dist = abs(z_i - z_j)
            if dist < min_dist:
                min_dist = dist
    return min_dist


