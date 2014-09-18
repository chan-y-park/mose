import cmath
from numerics import *
from kswcf import progeny_2
from parameters import *

class Trajectory:
    """
    The trajectory class.

    Attributes: coordinates, periods, degeneracy, phase, charge, parents, \
    boundary_condition, count (shared).
    Methods: evolve, terminate (?).
    Arguments of the instance: (initial_charge, degeneracy, phase, parents, \
    boundary_condition)
    """

    count = 0

    def __init__(self, initial_charge, degeneracy, phase, parents, 
                 boundary_condition, color='b'):
        self.degeneracy = degeneracy
        self.phase = phase
        self.parents = parents
        self.boundary_condition = boundary_condition
        self.initial_charge = initial_charge        # the initial charge
        self.color = color
        self.initial_point = None       # to be defined below
        self.evolve()
        Trajectory.count += 1 

    def __str__(self):
        return 'Trajectory info: initial charge %s , degeneracy %d, etc... ' % \
        (self.initial_charge, self.degeneracy)

    def get_color(self):
        return self.color
    
    def get_xcoordinates(self):
        return [z[0] for z in self.coordinates]

    def get_ycoordinates(self):
        return [z[1] for z in self.coordinates]

    def evolve(self):
        if verb: 
            print "\nEvolving trajectory %d " % Trajectory.count
        from numpy import concatenate
        if self.parents[0].__class__.__name__ == 'BranchPoint':
            ## for a *primary* kwall the boundary conditions are a bit particular:
            self.initial_point = self.parents[0]
            u0 = self.boundary_condition[0]
            sign = self.boundary_condition[1]
            pw_data = grow_primary_kwall(u0, sign, g2, g3, self.phase, primary_options)
            self.coordinates = pw_data[0]
            self.periods = pw_data[1]
            # now switch to picard-fuchs evolution
            bc = set_primary_bc(self)
            pw_data_pf = grow_pf(bc, g2, g3, self.phase, options)
            self.coordinates = concatenate(( self.coordinates, [ [row[0], row[1]] for row in pw_data_pf ] ))
            self.periods = concatenate(( self.periods , [row[2] + 1j* row[3] for row in pw_data_pf] ))
            self.check_cuts()
        elif self.parents[0].__class__.__name__ == 'Trajectory':
            ### recall that self.boundary_conditions in this case are set by set_bc(...)
            ### and they are formatted as [u0, eta0, d_eta0, intersection]
            ### the latter being an IntersectionPoint object
            self.initial_point = self.boundary_condition[3]
            pw_data_pf = grow_pf(self.boundary_condition[0:3], g2, g3, self.phase, options)
            self.coordinates =  [ [row[0], row[1]] for row in pw_data_pf ]
            self.periods =  [ row[2] + 1j* row[3] for row in pw_data_pf ] 
            self.check_cuts()

    def charge(self, point):
        return self.initial_charge
        # to be updated by taking into account branch-cuts, 
        # will use data in self.splittings

    def check_cuts(self):
        #print "Checking intersections with cuts, determining charges"
        self.splittings = (55, 107, 231) 
        self.local_charge = (self.initial_charge, (2,1), (0,-1), (1,1))
        # determine at which points the wall crosses a cut, for instance
        # (55,107,231) would mean that we change charge 3 times
        # hence self.splittings would have length, 3 while
        # self.local_charge would have length 4.
        # local charges are determined one the branch-cut data is given,
        # perhaps computed by an external function.



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
    """The BranchCut class.

    Attributes: locus, charge
    Arguments: branch-point (as an object), direction (as a phase e^(i phi))
    """
    from parameters import branch_cut_cutoff

    count = 0

    def __init__(self, branch_point, phase):
        self.charge = branch_point.charge
        self.locus = (branch_point.locus, 
                      complex(branch_point.locus + branch_cut_cutoff * phase))
        BranchCut.count += 1

    def __str__(self):
        return 'Branch cut info: charge %s, start-end-points %s ' % \
        (self.charge, self.locus)




class IntersectionPoint:
    """The IntersectionPoint class.

    Attributes: 
    locus (point on moduli space), 
    index_1 (index of intersection point within parent number 1, important to determine the charge at intersection), 
    index_2 (index of intersection point within parent number 2, important to determine the charge at intersection), 
    genealogy

    Arguments: 
    data (as a triplet of [u, index_1, index_2]), 
    parents (as list of trajectories, ie objects)    
    """

    def __init__(self,data,parents):

        self.parents = parents
        self.locus = data[0]
        self.index_1 = data[1]
        self.index_2 = data[2]
        self.charges = [parents[0].charge(self.index_1), parents[1].charge(self.index_2)]
        self.degeneracies = [parents[0].degeneracy, parents[1].degeneracy]
        self.genealogy = build_genealogy_tree(self)
        self.phase = parents[0].phase


class MarginalStabilityWall:
    """
    Th MS wall class.
    Attributes:
    charges, degeneracy -- those of (any) two K-walls meeting theta_range
    points -- the list of actual intersectionpoint objects
    locus -- the list of points (UNSORTED FOR NOW)
    """

    count = 0

    def __init__(self, all_intersections):
        self.charges = all_intersections[0].charges
        self.degeneracies = all_intersections[0].degeneracies
        self.points = all_intersections
        self.locus = [intersection.locus for intersection in all_intersections]
        MarginalStabilityWall.count += 1



def prepare_branch_locus(g2, g3, phase):
    """Find branch points and build branch cuts."""
    fixed_charges = [ [1, 0], [-1, 2] ] ### Must update with actual charge at branch-point
    branch_point_loci = map(complex, find_singularities(g2, g3))
    bpts = [BranchPoint(branch_point_loci[i], fixed_charges[i])  for i in range(len(branch_point_loci))]
    bcts = [BranchCut(bpts[i], phase) for i in range(len(bpts))]
    return [bpts, bcts]




def build_first_generation(branch_points, phase, g2, g3):
    """Construct the primary Kwalls"""
    traj = []
    Trajectory.count = 0
    bp = branch_points[0]
    colors = ['r','b','g','k']
    
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, +1],colors[i]))
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, -1],colors[i]))
    
    return traj




def new_intersections(kwalls,new_kwalls):
    """Find new wall-wall intersections"""

    from parameters import dsz_matrix, verb
    new_ints = []

    from parameters import intersection_range ## parameter related only to the temporary intersection algorithm
    
    if verb:
        print "\nsearching intersections of %s overall kwalls with %s new_kwalls" % (len(kwalls)+len(new_kwalls), len(new_kwalls))

    ### note: I am excluding some cases from being checked for intersections, see the if statements below
    for i_1, traj_1 in list(enumerate(kwalls)):
        for i_2, traj_2 in list(enumerate(new_kwalls)):
            if (dsz_pairing(traj_1.charge(0), traj_2.charge(0), dsz_matrix) != 0 and traj_1.parents != traj_2.parents and (not(traj_1 in traj_2.parents)) ):
                intersections_list = find_intersections(traj_1, traj_2)
                new_ints += [ IntersectionPoint(intersection, [traj_1, traj_2]) for intersection in intersections_list] 


    for i_1, traj_1 in list(enumerate(new_kwalls)):
        for i_2, traj_2 in list(enumerate(new_kwalls))[i_1+1 : ]:
            if (dsz_pairing(traj_1.charge(0), traj_2.charge(0), dsz_matrix) != 0 and traj_1.parents != traj_2.parents and not(traj_1 in traj_2.parents) and not(traj_2 in traj_1.parents)):
                intersections_list = find_intersections(traj_1, traj_2)
                new_ints += [ IntersectionPoint(intersection, [traj_1, traj_2]) for intersection in intersections_list] 

    if verb: 
        print "\nEvaluating intersections of NEW Kwalls with ALL Kwalls: found %d of them" % len(new_ints)
    return new_ints
    



def iterate(n,kwalls,new_kwalls,intersections):
    """Iteration"""
    from parameters import verb
    for i in range(n):
        new_ints = new_intersections(kwalls, new_kwalls)
        intersections += new_ints
        kwalls += new_kwalls
        if verb:
            print "\ncreating new trajectories"
        new_kwalls = build_new_walls(new_ints)

    return kwalls, new_kwalls, intersections




def build_new_walls(new_intersections):     
    """Build K-walls from new intersections"""
    new_walls = []
    
    from numpy import array

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
            charge = sibling[0] # this is the charge formatted wrt the basis of parent charges
            actual_charge = list( charge[0] * array(gamma_1) + charge[1] * array(gamma_2) )
            degeneracy = sibling[1]
            boundary_condition = set_bc(intersection, charge)
            new_walls.append( Trajectory(actual_charge, degeneracy, phase, parents, boundary_condition) )

    return new_walls
    



def build_genealogy_tree(intersection):
    """
    this function will return the genealogy tree of an intersection \
    in the form of a list such as [[BP1,BP2],BP3] for an intersection with the\
    obvious history
    """
    
    parents = intersection.parents
    index_1 = intersection.index_1
    index_2 = intersection.index_2
    
    ### determine who's mom and who's dad by relative orientation
    delta_z_1 = complexify(parents[0].coordinates[index_1+1]) - complexify(parents[0].coordinates[index_1])
    delta_z_2 = complexify(parents[1].coordinates[index_2+1]) - complexify(parents[1].coordinates[index_2])

    if cmath.phase(delta_z_1 / delta_z_2) > 0:
        dad = parents[0].initial_point
        mom = parents[1].initial_point
    else:
        dad = parents[1].initial_point
        mom = parents[0].initial_point

    return [dad.genealogy, mom.genealogy] 



def build_ms_walls(all_intersections):
    return "TO DO"
    ####### build MS walls based on 1) genealogies 2) charges of walls at the intersections



def phase_scan(theta_range):
    """Scans various values of theta, returns an array of IntersectionPoint objects.
    The argument is of the form: theta_range = [theta_in, theta_fin, steps]"""
    
    from parameters import g2, g3, primary_options, options, kwalls, new_kwalls, intersections, theta_cuts

    
    theta_in = theta_range[0]
    theta_fin = theta_range[1]
    steps = theta_range[2]

    angles = [ theta_in + i * (theta_fin - theta_in) / steps  for i in range(steps+1)]
    all_intersections = []
    all_kwalls = []

    branch_locus = prepare_branch_locus(g2, g3, theta_cuts)
    branch_points = branch_locus[0]
    branch_cuts = branch_locus[1]

    iter_count = 1

    for phase in angles:
        print "\n----------------------------------------------------------\
        \nIteration number %s: computing phase %s\
        \n----------------------------------------------------------" % (iter_count, phase)

        iter_count += 1
        new_kwalls = []
        kwalls = []
        intersections = []
        new_kwalls = build_first_generation(branch_points, phase, g2, g3)
        kwalls, new_kwalls, intersections = iterate(n_iter, kwalls, new_kwalls, 
            intersections)
        all_intersections += intersections
        all_kwalls += (kwalls + new_kwalls)

    ms_walls = build_ms_walls(all_intersections)

    return [all_intersections, all_kwalls, ms_walls]
