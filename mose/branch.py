class BranchPoint:
    """The BranchPoint class.

    Attributes: locus, charge
    Arguments: locus, charge
    """

    count = 0

    def __init__(self, locus, charge, monodromy_matrix, identifier):
        self.charge = charge
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
