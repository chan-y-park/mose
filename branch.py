

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
        BranchPoint.count += 1

    def __str__(self):
        return 'Branch point info: charge %s, locus %s ' % \
            (self.charge, self.locus)


class BranchCut:
    """
    The BranchCut class.

    Attributes: locus, charge
    Arguments: branch-point (as an object), direction (as a phase e^(i phi))
    """
    count = 0

    def __init__(self, branch_point, phase, branch_cut_cutoff):
        self.charge = branch_point.charge
        self.locus = (
            branch_point.locus,
            complex(branch_point.locus + branch_cut_cutoff * phase)
        )
        BranchCut.count += 1

    def __str__(self):
        return 'Branch cut info: charge %s, start-end-points %s ' % \
            (self.charge, self.locus)
