import numpy as np

class BranchPoint:
    """The BranchPoint class.

    Attributes: locus, charge
    Arguments: locus, charge
    """

    count = 0

    def __init__(self, locus, charge):
        self.charge = charge
        self.locus = locus
        self.count = BranchPoint.count
        self.genealogy = self
        # NEED TO IMPORT THE MONODROMY FROM WEIERSTRASS
        self.monodromy_matrix = np.identity(len(charge))
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
