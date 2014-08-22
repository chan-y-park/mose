import cmath
from numerics import *

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
                 boundary_condition):
        self.degeneracy = degeneracy
        self.phase = phase
        self.parents = parents
        self.boundary_condition = boundary_condition
        self.initial_charge = initial_charge        # the initial charge
        Trajectory.count += 1 
        # Make trajectory evolve automatically at this point?
        self.evolve()

    def __str__(self):
        return 'Trajectory info: initial charge %s , degeneracy %d, etc... ' % \
        (self.initial_charge, self.degeneracy)

    def get_xcoordinates(self):
        return [z[0] for z in self.coordinates]

    def get_ycoordinates(self):
        return [z[1] for z in self.coordinates]

    def evolve(self):
        if self.parents[0].__class__.__name__ == 'BranchPoint':
            u0 = self.boundary_condition[0]
            sign = self.boundary_condition[1]
            g2 = self.boundary_condition[2][0]
            g3 = self.boundary_condition[2][1]
            options = self.boundary_condition[2][2]
            self.coordinates = grow_primary_kwall(u0, sign, g2, g3, self.phase, options)

            #######     TO DO FROM HERE ON
            self.periods = (0, 2, 4)        # periods of the holomorphic one-form
            self.check_cuts()
        elif self.parents[0].__class__.__name__ == 'Trajectory':
            YYYY
            #print "Evolving trajectory"
            self.coordinates = (0+0j, 1+1j, 2+2j)       # points of the trajectory, 
            # must write algorithm to evolve
            self.periods = (0, 2, 4)        # periods of the holomorphic one-form
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
    cutoff = 10         # how far way from the singularity the locus of the 
                        # branch cut extends

    def __init__(self, branch_point, phase):
        self.charge = branch_point.charge
        self.locus = (branch_point.locus, 
                      branch_point.locus + BranchCut.cutoff * phase)
        BranchCut.count += 1

    def __str__(self):
        return 'Branch cut info: charge %s, start-end-points %s ' % \
        (self.charge, self.locus)




class IntersectionPoint:
    """The IntersectionPoint class.

    Attributes: locus, genealogy
    Arguments: locus, parents (as list of trajectories, ie objects), genealogy
    """

    def __init__(self,locus,parents):
        self.parents = parents
        self.locus = locus
        self.genealogy = build_genealogy_tree(parents)

    def __str__(self):
        return 'Intersection info: locus %s, parents (%s,%s) ' % \
        (self.locus, self.parents[0], self.parents[1])





def prepare_branch_locus(g2, g3, phase):
    """Find branch points and build branch cuts."""
    fixed_charge = (0,0) ### Must update with actual charge at branch-point
    branch_point_loci = map(complex, find_singularities(g2, g3))
    bpts = [BranchPoint(branch_point_loci[i], fixed_charge)  for i in range(len(branch_point_loci))]
    bcts = [BranchCut(bpts[i], phase) for i in range(len(bpts))]
    return [bpts, bcts]




def build_first_generation(branch_points, phase, g2, g3, options):
    """Construct the primary Kwalls"""
    #global new_kwalls
    #global kwalls
    #global intersections
    #new_kwalls = []
    #kwalls = []
    #intersections = []
    extra_parameters = [g2, g3, options]
    traj = []
    bp = branch_points[0]
    # here insert creation algorithm for the primary k-walls
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, +1, extra_parameters]))
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, -1, extra_parameters]))
    #[ Trajectory(bp.charge, 1, phase, bp, None) for i in range(len(branch_points))]
    #t1 = Trajectory(bp.charge, 1, phase, bp, None)
    #t2 = Trajectory((-1, 2), 1, 0, (23, 7), "bc1")
    #
    return traj
    #kwalls = [t1, t2]
    #new_kwalls = [t1, t2]
    #print "\nConstructed primary Kwalls, there are %d of them" % len(new_kwalls)





def new_intersections():
    """Find new wall-wall intersections"""
    # here insert algorithm that computes all intersections of NEW kwalls with 
    # ALL kwalls
    int1 = IntersectionPoint(0+0j, [new_kwalls[0], new_kwalls[1]])
    new_ints = [int1]
    print "Evaluating intersections of NEW Kwalls with ALL Kwalls: found %d of them" % len(new_ints)
    return new_ints
    



def iterate(n):
    global intersections
    global kwalls
    global new_kwalls
    """Iteration"""
    new_ints = new_intersections()
    new_kwalls = build_new_walls(new_ints)
    kwalls = kwalls + new_kwalls
    intersections = intersections + new_ints    
    # Make it repeat n times, once it works




def build_new_walls(intersections):     
    #the argument here will be provided directly by function new_intersections()
    """Build K-walls from new intersections"""
    global t3
    print "Constructing new walls"
    # here insert algorithm that computes the new walls
    t3 = Trajectory((0,2), 1, 0, (12,35), "bc3")
    new_walls = [t3]
    #
    return new_walls
    



def build_genealogy_tree(parents):
    return "this function will return the genealogy tree of an intersection \
    in the form of a list such as [[BP1,BP2],BP3] for an intersection with the\
    obvious history"
    # Genealogy is crucial for identifying which points should belong to the 
    # same MS wall.
    # It will not be necessary to develop this until everything else is in place.






###### Some Checks ######


### An Example of a Workflow


#print "Preparing the branch locus:"
#phase = cmath.exp(1.23j)
#prepare_branch_locus(phase)

#print "\nPreparing the primary kwalls:"

#build_first_generation()

#print "\n the kwalls are:"
#print kwalls[0]
#print kwalls[1]

#print "\n the new kwalls are:"
#print new_kwalls[0]
#print new_kwalls[1]

#print "\n the intersections are:"
#print intersections

#print "\n --going to the next generation-- \n"
#iterate(1)

#print "\n the kwalls are:"
#print kwalls[0]
#print kwalls[1]
#print kwalls[2]

#print "\n the new kwalls are:"
#print new_kwalls[0]

#print "\n the intersections are:"
#print intersections[0]


### Some Command Tests

#print t1
#print t2


#print "coordinates: %s" % (t1.coordinates,)
#print "periods: %s" % (t1.periods,)
#print "charge at step %d: %s" % (55,t1.charge(55))



#print BC1
#print BP1

#print "There are %d trajectories." % Trajectory.count
#print "There are %d branch points." % BranchPoint.count
#print "There are %d branch cuts." % BranchCut.count

#print t1.__class__.__name__
#print BP1.__class__.__name__
#print BC1.__class__.__name__



#print intersections[0].genealogy
