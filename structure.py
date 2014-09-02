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
            g2 = self.boundary_condition[2][0]
            g3 = self.boundary_condition[2][1]
            primary_options = self.boundary_condition[2][2]
            options = self.boundary_condition[2][3]
            pw_data = grow_primary_kwall(u0, sign, g2, g3, self.phase, primary_options)
            self.coordinates = pw_data[0]
            self.periods = pw_data[1]
            # 
            bc = set_primary_bc(self)
            pw_data_pf = grow_pf(bc, g2, g3, self.phase, options)
            self.coordinates = concatenate(( self.coordinates, [ [row[0], row[1]] for row in pw_data_pf ] ))
            self.periods = concatenate(( self.periods , [row[2] + 1j* row[3] for row in pw_data_pf] ))
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




def build_first_generation(branch_points, phase, g2, g3, primary_options, options):
    """Construct the primary Kwalls"""
    
    extra_parameters = [g2, g3, primary_options, options]
    traj = []
    Trajectory.count = 0
    bp = branch_points[0]
    colors = ['r','b','g','k']
    
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, +1, extra_parameters],colors[i]))
    for i in range(len(branch_points)):
        bp = branch_points[i]
        traj.append(Trajectory(bp.charge, 1, phase, [bp], [bp.locus, -1, extra_parameters],colors[i+2]))
    
    return traj




def new_intersections(kwalls,new_kwalls):
    """Find new wall-wall intersections"""
    # here insert algorithm that computes all intersections of NEW kwalls with 
    # ALL kwalls
    int1 = IntersectionPoint(0+0j, [new_kwalls[0], new_kwalls[1]])
    new_ints = [int1]
    print "Evaluating intersections of NEW Kwalls with ALL Kwalls: found %d of them" % len(new_ints)
    return new_ints
    



def iterate(n,kwalls,new_kwalls,intersections):
    """Iteration"""
    new_ints = new_intersections(kwalls, new_kwalls)
    intersections = intersections + new_ints    
    new_kwalls = build_new_walls(new_ints)
    kwalls = kwalls + new_kwalls
    # Make it repeat n times, once it works




def build_new_walls(intersections):     
    #the argument here will be provided directly by function new_intersections()
    """Build K-walls from new intersections"""
    # global t3
    print "Constructing new walls"
    # here insert algorithm that computes the new walls: 
    # includes the KSWCF to see what new kwalls to create, 
    # and must receive the appropriate boundary conditions for PF evolution
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

