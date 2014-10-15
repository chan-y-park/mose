import logging
import multiprocessing
from numpy import array
from k_wall import PrimaryKWall, DescendantKWall
from intersection import HitTable
from intersection_point import find_new_intersections
from misc import complexify
from kswcf import progeny_2
from k_wall import KWall

class KWallNetwork:
    def __init__(self, theta, fibration, bin_size):
        self.phase = theta
        self.fibration = fibration
        self.hit_table = HitTable(bin_size)
        self.k_walls = []
        self.intersections = []
        KWall.count = 0

    def grow(self, primary_nint_range, nint_range,
             trajectory_singularity_threshold, pf_odeint_mxstep, 
             n_iterations, ks_filtration_degree):
        logging.info('Growing a K-Wall network at phase %.8f', self.phase)
        ##############################
        # First, grow primary k-walls.
        ##############################
        primary_k_walls = []

        for sign in [+1, -1]:
            for bp in self.fibration.branch_points:
                # logging.info('Evolving primary K-wall #%d', 
                                # len(primary_k_walls))
                k_wall = PrimaryKWall(
                    list(sign * array(bp.charge)),      #initial_charge 
                    1,              #degeneracy 
                    self.phase, 
                    [bp],           #parents 
                    self.fibration,
                    [bp.locus, sign], #boundary_condition
                    primary_nint_range
                )
                k_wall.evolve(nint_range, trajectory_singularity_threshold,
                                pf_odeint_mxstep)
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
            # print "flag 1 %s" % len(primary_k_walls)
            logging.info('Iteration #%d', i+1 )
            logging.debug('len(k_walls) = %d', len(self.k_walls))
            new_intersections = find_new_intersections(
                self.k_walls, new_k_walls, self.intersections, 
                self.hit_table, self.fibration.dsz_matrix
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
                    self.fibration.dsz_matrix, 
                    ks_filtration_degree
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
                    k_wall.evolve(nint_range, trajectory_singularity_threshold,
                                    pf_odeint_mxstep)
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

def parallel_grow(k_wall_network, *args):
    k_wall_network.grow(*args)

def construct_k_wall_networks(fibration, bin_size, primary_nint_range,
                              nint_range, trajectory_singularity_threshold,
                              pf_odeint_mxstep, n_iterations,
                              ks_filtration_degree, theta_range, n_process=0):
    """
    Scans various values of theta, returns an array of 
    IntersectionPoint objects.
    The argument is of the form: theta_range = [theta_i, theta_f, steps]
    """
    
    theta_i, theta_f, steps = theta_range

    angles = [theta_i + i*(theta_f - theta_i)/steps for i in range(steps)]
    k_wall_networks = []

    for phase in angles:
        kwn = KWallNetwork(phase, fibration, bin_size)
        k_wall_networks.append(kwn)

    manager = multiprocessing.Manager()
    shared_kwns = manager.list(k_wall_networks)

#    for kwn in shared_kwns:
#        p = multiprocessing.Process(
#                target=kwn.grow,
#                args=(primary_nint_range, nint_range,
#                      trajectory_singularity_threshold, pf_odeint_mxstep,
#                      n_iterations, ks_filtration_degree)
#        )
#        p.start()
#        p.join()

    if(n_process == 0):
        n_process = multiprocessing.cpu_count()
    logging.debug('Number of processes: %d', n_process)
    pool = multiprocessing.Pool(processes=n_process)
    results = [pool.apply_async(parallel_grow,
                                args=(kwn, primary_nint_range, nint_range,
                                      trajectory_singularity_threshold,
                                      pf_odeint_mxstep, n_iterations,
                                      ks_filtration_degree))
               for kwn in shared_kwns]
    pool.close()
    pool.join()
#    for p in results:
#        p.wait()
    return k_wall_networks
