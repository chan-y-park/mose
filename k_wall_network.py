from numpy import array
from k_wall import KWall
from intersection import HitTable
from misc import complexify
from kswcf import progeny_2

class KWallNetwork:
    def __init__(self, theta, fibration, max_range, bin_size):
        self.phase = theta
        self.fibration = fibration
        self.hit_table = HitTable(max_range, bin_size)
        self.k_walls = []
        self.intersections = []

    def grow(self, n_iterations):
        ##############################
        # First, grow primary k-walls.
        ##############################
        primary_k_walls = []

        for bp in self.fibration.branch_points:
            k_wall = KWall(initial_charge = bp.charge, degeneracy = 1, 
                            phase = self.phase, parents = [bp], 
                            boundary_condition = [bp.locus, +1])
            if (not k_wall.singular):
                primary_k_walls.append(k_wall)
            else:
                logging.info(
                                '**************\n'
                                'SINGULAR K-WALL! WILL BE DROPPED\n'
                                '**************\n'
                            )

            k_wall = KWall(initial_charge = bp.charge, degeneracy = 1, 
                            phase = self.phase, parents = [bp], 
                            boundary_condition = [bp.locus, -1])
            if (not k_wall.singular):
                primary_k_walls.append(k_wall)
            else:
                logging.info(
                                '**************\n'
                                'SINGULAR K-WALL! WILL BE DROPPED\n'
                                '**************\n'
                            )

        #############################
        # Now get descendant k-walls.
        #############################

        new_k_walls = primary_k_walls 
        for i in range(n_iterations):
            new_intersections = find_new_intersections(
                self.k_walls, new_k_walls, self.intersections, self.hit_table
            )
            self.intersections += new_ints
            self.k_walls += new_k_walls
            logging.info('creating new K-walls from intersections.')
            new_k_walls = build_k_walls_from_intersections(new_intersections)

def build_k_walls_from_intersections(new_intersections):     
    """Build K-walls from new intersections"""
    new_k_walls = []
    
    for intersection in new_intersections:
        parents = intersection.parents
        gamma_1 = parents[0].charge(intersection.index_1)
        gamma_2 = parents[1].charge(intersection.index_2)
        omega_1 = parents[0].degeneracy
        omega_2 = parents[1].degeneracy
        u_0 = intersection.locus
        phase = intersection.phase
        progeny = progeny_2([[gamma_1, omega_1], [gamma_2, omega_2]])
        for sibling in progeny:
            # the charge formatted wrt the basis of parent charges
            charge, degeneracy = sibling 
            actual_charge = list(charge[0]*array(gamma_1) + 
                                    charge[1]*array(gamma_2))
            boundary_condition = set_boundary_condition(intersection, charge)
            k_wall = KWall(actual_charge, degeneracy, phase, 
                                parents, boundary_condition)
            if (not k_wall.singular):
                new_k_walls.append(k_wall)
            else:
                logging.info(
                                '**************\n'
                                'SINGULAR K-WALL! WILL BE DROPPED\n'
                                '**************\n'
                            )
    return new_k_walls
    
def set_boundary_condition(intersection, charge):
    """ 
    Employs data of class KWall to produce correctly formatted 
    boundary conditions for grow_pf().

    It also passes the argument "intersection" through, this will be employed
    in the KWall class to keep track of genealogical data.

    Arguments: (intersection, charge)
        intersection: must be an instance of the IntersecionPoint class.
        charge: must be the charge relative to the parents' charge basis, 
            hence a list of length 2.
    """
    u_0 = intersection.locus

    parents = intersection.parents
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    path_1 = map(complexify, parents[0].coordinates)
    path_2 = map(complexify, parents[1].coordinates)
    periods_1 = parents[0].periods
    periods_2 = parents[1].periods
    eta_1 = periods_1[index_1]
    eta_2 = periods_2[index_2]
    d_eta_1 = ((periods_1[index_1] - periods_1[index_1-1]) / 
                (path_1[index_1] - path_1[index_1-1]))
    d_eta_2 = ((periods_2[index_2] - periods_2[index_2-1]) / 
                (path_2[index_2] - path_2[index_2-1]))
    # NOTE: The use of complex() is necessary here, because sometimes 
    # the charge vector wil be deriving from an algorithm using sympy, 
    # and will turn j's into I's...
    eta_0 = eta_1*complex(charge[0]) + eta_2*complex(charge[1])  
    d_eta_0 = d_eta_1*complex(charge[0]) + d_eta_2*complex(charge[1])

    return [u_0, eta_0, d_eta_0, intersection]

def construct_k_wall_networks(fibration, theta_range, n_iter):
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
        Iteration number {}: computing phase {}
        ----------------------------------------------------------
        """.format(iter_count, phase)

        kwn = KWallNetwork(phase, fibration)

        kwn.grow(n_iter)
        k_wall_networks.append(kwn)
        # kwallplot(kwalls+new_kwalls, INTERSECTION_SEARCH_RANGE,
        #             branch_points=branch_points) 
        #kwallplot(kwalls+new_kwalls, INTERSECTION_SEARCH_RANGE,
        #            bin_size=INTERSECTION_SEARCH_BIN_SIZE,
        #            intersection_points=intersections,
        #            hit_table=ht, mark_data_plot=True) 


    return k_wall_networks
