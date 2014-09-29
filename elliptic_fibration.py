import logging
import sympy as sym
from branch import BranchPoint, BranchCut

class EllipticFibration:
    def __init__(self, g2, g3, branch_point_charges,
                    branch_cut_phase, branch_cut_cutoff):

        self.g2 = g2
        self.g3 = g3
        branch_point_loci = map(complex, find_singularities(g2, g3))
        self.branch_points = [
            BranchPoint(branch_point_loci[i], branch_point_charges[i])  
            for i in range(len(branch_point_loci))
        ]
        self.branch_cuts = [
            BranchCut(bp, branch_cut_phase, branch_cut_cutoff) 
            for bp in self.branch_points
        ]

def find_singularities(g2, g3):
    """
    find the singularities on the Coulomb branch
    """
    u = sym.Symbol('u')
    discriminant = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
    logging.info('discriminant: %s', discriminant)
    logging.info('singularities: %s', sym.solve(discriminant, u))
    return sym.solve(discriminant, u)


