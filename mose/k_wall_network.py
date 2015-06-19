import logging
import signal
import multiprocessing
# import pdb
from multiprocessing.managers import BaseManager
from numpy import array
from elliptic_fibration import EllipticFibration
from k_wall import PrimaryKWall, DescendantKWall
from hit_table import HitTable
from intersection_point import find_new_intersection_points
from intersection_point import find_new_intersections_using_hit_table
from misc import complexify, cut_singular_k_wall, id_generator
from kswcf import progeny_2
from k_wall import KWall



class KWallNetwork:
    def __init__(self, phase, fibration, config, use_hit_table=False):
        self.phase = phase
        self.fibration = fibration
        if use_hit_table is True:
            self.hit_table = HitTable(
                config['intersection search']['bin_size']
            )
        else:
            self.hit_table = None
        self.k_walls = []
        self.intersections = []


    def save_data(self, file_object, **kwargs):
        pass


    def grow(self, config):
        #primary_nint_range = config['ODE']['primary_k_wall_ode_range']
        #nint_range = config['ODE']['ode_range']
        trajectory_singularity_threshold = (
            config['ODE']['trajectory_singularity_threshold']
        )
        ode_size_of_step = config['ODE']['size_of_step']
        ode_num_steps = config['ODE']['num_steps']
        n_iterations = config['K-wall network']['n_iterations']
        ks_filtration_degree = config['KSWCF']['filtration_degree']

        logging.info('Growing a K-Wall network at phase %.8f', self.phase)
        ##############################
        # First, grow primary k-walls.
        ##############################
        primary_k_walls = []

        for sign in [1.0, -1.0]:
            ### Here the sign list distinguishes between \gamma and -\gamma
            ### both of whose K-walls are generated from the same discriminant
            ### locus. The distinction between these is determined ultimately 
            ### by the Weierstrass analysis, within the class PrimaryKWall
            for bp in self.fibration.branch_points:
            #for idx, bp in enumerate(self.fibration.branch_points):
                logging.info(
                    'Evolving primary K-wall #%d', len(primary_k_walls)
                )
                # NOTE: the assignment of charge is only provisional,
                # within the creation of the Kwall we will determine the 
                # charge from the parent branch point, taking care of 
                # the distinction between \gamma and -\gamma.
                k_wall = PrimaryKWall(
                    degeneracy=1,
                    phase=self.phase,
                    parents=[bp],
                    fibration=self.fibration,
                    initial_condition=[bp.locus, sign],
                    #label="K-wall #{}".format(len(self.k_walls))
                    identifier=id_generator(),
                    network=self
                )
                k_wall.grow_pf(
                    trajectory_singularity_threshold,
                    ode_size_of_step,   
                    ode_num_steps,
                )
                if (not k_wall.singular):
                    primary_k_walls.append(k_wall)
                else:
                    ### NOT CUTTING ANY MORE
                    # cut_singular_k_wall(k_wall)
                    primary_k_walls.append(k_wall)
                    # logging.info(
                    #     """
                    #     **************
                    #     SINGULAR K-WALL! WILL BE DROPPED.
                    #     **************
                    #     """
                    # )

        #############################
        # Now grow descendant k-walls.
        #############################

        new_k_walls = primary_k_walls
        for i in range(n_iterations):
            logging.info('Iteration #%d', i+1)
            logging.debug('len(k_walls) = %d', len(self.k_walls))
            if self.hit_table is None:
                new_intersections = find_new_intersection_points(
                    self.k_walls, new_k_walls, self.intersections,
                    self.fibration.dsz_matrix,
                )
            else:
                new_intersections = find_new_intersections_using_hit_table(
                    self.k_walls, new_k_walls, self.intersections,
                    self.hit_table, self.fibration.dsz_matrix,
                )
            self.intersections += new_intersections
            self.k_walls += new_k_walls
            new_k_walls = []
            logging.info('Creating new K-walls from intersections.')

            # Build K-walls from new intersections.

            for intersection in new_intersections:
                ### The parents are sorted by the IntersectionPoint class
                ### in such a way that the KSWCF is K_2 K_1 = K_1 ... K_2
                ### i.e. Arg(Z_1) < Arg(Z_2) on the LHS
                parents = intersection.parents
                gamma_1 = parents[0].charge(intersection.indices[0])
                gamma_2 = parents[1].charge(intersection.indices[1])
                omega_1 = parents[0].degeneracy
                omega_2 = parents[1].degeneracy
                u_0 = intersection.locus
                phase = intersection.phase
                progeny = progeny_2(
                    [[gamma_1, omega_1], [gamma_2, omega_2]],
                    self.fibration.dsz_matrix,
                    ks_filtration_degree,
                )
                logging.debug('progeny = %s', progeny)
                if progeny == None:
                    pass
                else:
                    for sibling in progeny:
                        # the charge formatted wrt the basis of parent charges
                        charge, degeneracy = sibling
                        actual_charge = list(
                            charge[0]*array(gamma_1) + charge[1]*array(gamma_2)
                        )
                        k_wall = DescendantKWall(
                            initial_charge=actual_charge,
                            degeneracy=degeneracy,
                            phase=phase,
                            parents=parents,
                            fibration=self.fibration,
                            intersection=intersection,
                            charge_wrt_parents=charge,
                            #label="K-wall #{}".format(len(self.k_walls))
                            identifier=id_generator(),
                            network=self
                        )
                        k_wall.grow_pf(
                            trajectory_singularity_threshold,
                            ode_size_of_step,   
                            ode_num_steps,
                        )
                        if (not k_wall.singular):
                            new_k_walls.append(k_wall)
                        else:
                            ### No longer cutting kwalls!
                            ### we approximate the bad region with a straight
                            ### line now
                            # cut_singular_k_wall(k_wall)
                            new_k_walls.append(k_wall)

        ### End of iterations.
        self.k_walls += new_k_walls

class FibrationManager(BaseManager): 
    pass

FibrationManager.register('fibration', EllipticFibration)

def init_process():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def parallel_get_k_wall_network(phase, fibration, config, 
                                shared_n_started_kwns, shared_n_finished_kwns):
    shared_n_started_kwns.value += 1
    logging.info('Start generating K-wall network #%d',
                 shared_n_started_kwns.value)
    kwn = KWallNetwork(phase, fibration, config,)
    kwn.grow(config)
    shared_n_finished_kwns.value += 1
    logging.info('Finished generating K-wall network #%d',
                 shared_n_finished_kwns.value)
    return kwn


def construct_k_wall_networks(fibration, config):
    """
    Scans various values of theta, returns an array of
    IntersectionPoint objects.
    The argument is of the form: theta_range = [theta_i, theta_f, steps]
    """
 
    theta_range = config['MS wall']['theta_range']
    n_processes= config['multiprocessing']['n_processes']

    logging.debug('theta_range = %s', theta_range)
    theta_i, theta_f, steps = theta_range

    angles = [theta_i + i*(theta_f - theta_i)/steps for i in range(steps)]
    k_wall_networks = []

    manager = multiprocessing.Manager()
    shared_n_started_kwns = manager.Value('i', 0)
    shared_n_finished_kwns = manager.Value('i', 0)

#    for kwn in shared_kwns:
#        p = multiprocessing.Process(
#                target=kwn.grow,
#                args=(primary_nint_range, nint_range,
#                      trajectory_singularity_threshold, pf_odeint_mxstep,
#                      n_iterations, ks_filtration_degree)
#        )
#        p.start()
#        p.join()

    logging.debug('n_processes: %d', n_processes)
    if(n_processes == 0):
        n_processes = multiprocessing.cpu_count()
    logging.debug('Number of processes: %d', n_processes)
    pool = multiprocessing.Pool(n_processes, init_process)
    try:
        results = [
            pool.apply_async(
                parallel_get_k_wall_network,
                args=(phase, fibration, config, shared_n_started_kwns,
                      shared_n_finished_kwns,)
            ) for phase in angles
        ]
        pool.close()

        for r in results:
            k_wall_networks.append(r.get())
            logging.info('job progress: %d/%d finished.',
                         shared_n_finished_kwns.value, len(angles))
    except KeyboardInterrupt:
        logging.warning('Caught ^C, terminating processes.')
        pool.terminate()
        pool.join()

    return k_wall_networks
