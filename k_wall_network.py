import logging
from numpy import array
from k_wall import PrimaryKWall, DescendantKWall
from intersection import HitTable
from intersection_point import find_new_intersections
from misc import complexify
from kswcf import progeny_2
from config import DSZ_MATRIX

class KWallNetwork:
    def __init__(self, theta, fibration, max_range, bin_size):
        self.phase = theta
        self.fibration = fibration
        self.hit_table = HitTable(max_range, bin_size)
        self.k_walls = []
        self.intersections = []

    def grow(self, primary_nint_range, nint_range, n_iterations):
        ##############################
        # First, grow primary k-walls.
        ##############################
        primary_k_walls = []

        for bp in self.fibration.branch_points:
            logging.info('Evolving primary K-wall #%d', 
                            len(primary_k_walls))
            k_wall = PrimaryKWall(
                bp.charge,      #initial_charge 
                1,              #degeneracy 
                self.phase, 
                [bp],           #parents 
                self.fibration,
                [bp.locus, +1], #boundary_condition
                primary_nint_range
            )
            k_wall.evolve(nint_range)
            if (not k_wall.singular):
                primary_k_walls.append(k_wall)
            else:
                logging.info(
                    """
                    **************
                    SINGULAR K-WALL! WILL BE DROPPED.
                    **************
                    """
                )

            logging.info('Evolving primary K-wall #%d', 
                            len(primary_k_walls))
            k_wall = PrimaryKWall(
                bp.charge,      #initial_charge 
                1,              #degeneracy 
                self.phase, 
                [bp],           #parents 
                self.fibration,
                [bp.locus, -1], #boundary_condition
                primary_nint_range
            )
            k_wall.evolve(nint_range)
            if (not k_wall.singular):
                primary_k_walls.append(k_wall)
            else:
                logging.info(
                    """
                    **************
                    SINGULAR K-WALL! WILL BE DROPPED.
                    **************
                    """
                )

        #############################
        # Now grow descendant k-walls.
        #############################

        new_k_walls = primary_k_walls 
        for i in range(n_iterations):
            logging.info('Iteration #%d', i)
            logging.debug('len(k_walls) = %d', len(self.k_walls))
            new_intersections = find_new_intersections(
                self.k_walls, new_k_walls, self.intersections, self.hit_table
            )
            self.intersections += new_intersections
            self.k_walls += new_k_walls
            new_k_walls = []
            logging.info('Creating new K-walls from intersections.')

            # Build K-walls from new intersections.
        
            for intersection in new_intersections:
                parents = intersection.parents
                gamma_1 = parents[0].charge(intersection.index_1)
                gamma_2 = parents[1].charge(intersection.index_2)
                omega_1 = parents[0].degeneracy
                omega_2 = parents[1].degeneracy
                u_0 = intersection.locus
                phase = intersection.phase
                progeny = progeny_2(
                    [[gamma_1, omega_1], [gamma_2, omega_2]],
                    DSZ_MATRIX
                )
                logging.debug('progeny = %s', progeny)
                for sibling in progeny:
                    # the charge formatted wrt the basis of parent charges
                    charge, degeneracy = sibling 
                    actual_charge = list(
                        charge[0]*array(gamma_1) + charge[1]*array(gamma_2)
                    )
                    k_wall = DescendantKWall(
                        actual_charge, #initial_charge
                        degeneracy, 
                        phase, 
                        parents, 
                        self.fibration,
                        intersection,
                        charge
                    )
                    k_wall.evolve(nint_range)
                    if (not k_wall.singular):
                        new_k_walls.append(k_wall)
                    else:
                        logging.info(
                            """
                            **************
                            SINGULAR K-WALL! WILL BE DROPPED.
                            **************
                            """
                        )
        # End of iterations.
        self.k_walls += new_k_walls

def construct_k_wall_networks(fibration, max_range, bin_size,
                                primary_nint_range, nint_range, n_iterations,
                                theta_range):
    """
    Scans various values of theta, returns an array of 
    IntersectionPoint objects.
    The argument is of the form: theta_range = [theta_i, theta_f, steps]
    """
    
    theta_i, theta_f, steps = theta_range

    angles = [theta_i + i*(theta_f - theta_i)/steps for i in range(steps)]
    k_wall_networks = []

    for phase in angles:
        print """
        ----------------------------------------------------------
        Network number {}: computing phase {}
        ----------------------------------------------------------
        """.format(len(k_wall_networks), phase)

        kwn = KWallNetwork(phase, fibration, max_range, bin_size)

        kwn.grow(primary_nint_range, nint_range, n_iterations)
        k_wall_networks.append(kwn)

    return k_wall_networks
