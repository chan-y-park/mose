import cmath
from numerics import *
from kswcf import *
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
        # Make trajectory evolve automatically at this point?
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
        print "\nEvolving trajectory %d " % Trajectory.count
        from numpy import concatenate
        if self.parents[0].__class__.__name__ == 'BranchPoint':
            ## for a *primary* kwall the boundary conditions are a bit particular:
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
            pw_data_pf = grow_pf(self.boundary_condition, g2, g3, self.phase, options)
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
        BranchPoint.count += 1 

    def __str__(self):
        return 'Branch point info: charge %s, locus %s ' % \
        (self.charge, self.locus)




class BranchCut:
    """The BranchCut class.

    Attributes: locus, charge
    Arguments: branch-point (as an object), direction (as a phase e^(i phi))
    """

    count = 0
    cutoff = 10.0         # how far way from the singularity the locus of the 
                        # branch cut extends

    def __init__(self, branch_point, phase):
        self.charge = branch_point.charge
        self.locus = (branch_point.locus, 
                      complex(branch_point.locus + BranchCut.cutoff * phase))
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
        self.genealogy = build_genealogy_tree(parents)

    def __str__(self):
        return 'Intersection info: locus %s, parents (%s,%s) ' % \
        (self.locus, self.parents[0], self.parents[1])



def prepare_branch_locus(g2, g3, phase):
    """Find branch points and build branch cuts."""
    fixed_charges = [ [1, 0], [-1, 2] ] ### Must update with actual charge at branch-point
    branch_point_loci = map(complex, find_singularities(g2, g3))
    bpts = [BranchPoint(branch_point_loci[i], fixed_charges[i])  for i in range(len(branch_point_loci))]
    bcts = [BranchCut(bpts[i], phase) for i in range(len(bpts))]
    return [bpts, bcts]




def build_first_generation(branch_points, phase, g2, g3, primary_options, options):
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

    from parameters import dsz_matrix
    new_ints = []

    from parameters import intersection_range ## parameter related only to the temporary intersection algorithm
    
    for i_1, traj_1 in list(enumerate(kwalls)):
        for i_2, traj_2 in list(enumerate(new_kwalls)):
            if (dsz_pairing(traj_1.charge(0), traj_2.charge(0), dsz_matrix) != 0 and traj_1.parents != traj_2.parents):
                intersections_list = find_intersections(traj_1, traj_2)
                # print "\nclosest point between trajectories [%s, %s]: %s" % (i_1, i_2, intersections_list)
                new_ints += [ IntersectionPoint(intersection, [traj_1, traj_2]) for intersection in intersections_list] 


    for i_1, traj_1 in list(enumerate(new_kwalls)):
        for i_2, traj_2 in list(enumerate(new_kwalls))[i_1+1 : ]:
            if (dsz_pairing(traj_1.charge(0), traj_2.charge(0), dsz_matrix) != 0 and traj_1.parents != traj_2.parents):
                intersections_list = find_intersections(traj_1, traj_2)
                # print "\nclosest point between trajectories [%s, %s]: %s" % (i_1, i_2, intersections_list)
                new_ints += [ IntersectionPoint(intersection, [traj_1, traj_2]) for intersection in intersections_list] 

    print "\nEvaluating intersections of NEW Kwalls with ALL Kwalls: found %d of them" % len(new_ints)
    return new_ints
    



def iterate(n,kwalls,new_kwalls,intersections):
    """Iteration"""
    for i in range(n):
        new_ints = new_intersections(kwalls, new_kwalls)
        intersections += new_ints
        kwalls += new_kwalls
        new_kwalls = build_new_walls(new_ints)

    return kwalls, new_kwalls, intersections




def build_new_walls(new_intersections):     
    """Build K-walls from new intersections"""
    
    print "\nConstructing new walls"
    new_walls = []
    
    from parameters import theta
    from numpy import array

    for intersection in new_intersections:
        parents = intersection.parents
        gamma_1 = parents[0].charge(intersection.index_1)
        gamma_2 = parents[1].charge(intersection.index_2)
        omega_1 = parents[0].degeneracy
        omega_2 = parents[1].degeneracy
        u_0 = intersection.locus
        # print "gamma_1 = %s" % gamma_1
        # print "gamma_2 = %s" % gamma_2
        # print "omega_1 = %s" % omega_1
        # print "omega_2 = %s" % omega_2
        # print "u_0 = %s" % u_0
   
        progeny = progeny_2([[gamma_1, omega_1], [gamma_2, omega_2]])
        for sibling in progeny:
            charge = sibling[0] # this is the charge formatted wrt the basis of parent charges
            actual_charge = list( charge[0] * array(gamma_1) + charge[1] * array(gamma_2) )
            degeneracy = sibling[1]
            boundary_condition = set_bc(intersection, charge)
            # print "actual_charge, , degeneracy, theta, parents, boundary_condition = %s " % [actual_charge, degeneracy, theta, parents, boundary_condition]
            new_walls.append( Trajectory(actual_charge, degeneracy, theta, parents, boundary_condition) )

    return new_walls
    



def build_genealogy_tree(parents):
    return "this function will return the genealogy tree of an intersection \
    in the form of a list such as [[BP1,BP2],BP3] for an intersection with the\
    obvious history"
    # Genealogy is crucial for identifying which points should belong to the 
    # same MS wall.
    # It will not be necessary to develop this until everything else is in place.

