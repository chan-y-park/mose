import sympy as sym

class EllipticFibration:
    def __init__(g2, g3, theta_cuts):
        self.branch_points, self.branch_cuts = \
            prepare_branch_locus(g2, g3, theta_cuts)


def prepare_branch_locus(g2, g3, phase):
    """
    Find branch points and build branch cuts.
    """
    # TODO: Must update with actual charge at branch-point
    fixed_charges = [[1, 0], [-1, 2]] 

    branch_point_loci = map(complex, find_singularities(g2, g3))
    bpts = [BranchPoint(branch_point_loci[i], fixed_charges[i])  
                for i in range(len(branch_point_loci))]
    bcts = [BranchCut(bpts[i], phase) for i in range(len(bpts))]

    return [bpts, bcts]
    
def find_singularities(g2, g3):
    """
    find the singularities on the Coulomb branch
    """

    u = sym.Symbol('u')
    discriminant = sym.simplify(g2 ** 3 - 27 * g3 ** 2)
    logging.info('discriminant: %s', discriminant)
    logging.info('singularities: %s', sym.solve(discriminant, u))
    return sym.solve(discriminant, u)


