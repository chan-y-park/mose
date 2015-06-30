"""
The main players of this package are the class

WeierstrassModelWithPaths

and the function

monodromy_at_point_wmodel



In the simplest use, a class for the Weierstrass model

y**2 = x**3 + f x + g

can be constructed by

WeierstrassModelWithPaths(fcoeff,gcoeff)

where fcoeff and gcoeff are arrays of coefficients of
the polynomials f and g, respectively. For example, when 

f = a x**2 + b x + c
g = d x + e

the arrays are given by

fcoeff = [a,b,c]
gcoeff = [d,e]

The object, wmodel, for this Weierstrass model can be
made by declaring

wmodel = WeierstrassModelWithPaths([a,b,c],[d,e])



The function monodromy_at_point_wmodel(n,w) returns the
monodromy matrix for the n'th discriminant locus of Weierstrass
model w. Here, when the Weierstrass model w has D discriminant
loci, n can be in the range [0,D-1]. The discriminant loci are
aligned in ascending order of the real part.

"""

import cmath
import logging
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from numpy import linalg as LA
from scipy.integrate import odeint
# from scipy.special import ellipk
from sympy.mpmath import ellipk
from itertools import combinations
from misc import period_A, period_B, real_part, rotate_poly, get_real_part



NEGLIGIBLE_BOUND = 0.1**4

U_PLANE_ROTATION_STEPS = 6

### Parameter for rotating the u-plane.
### Let d_min be the minimal distance |u_i - u_j|
### between discriminant loci.
### Then the u-plane will be rotated if the minimal
### separation between the REAL PARTS of any two loci
### is less than d_min * MINIMUM_DISCRIMINANT_LOCUS_SEPARATION
MINIMUM_DISCRIMINANT_LOCUS_SEPARATION = 0.1

#### Classes of Weierstrass Models ####

class WeierstrassModel:
    """
    WeierstrassModel(fcoeff,gcoeff)
    
    is a class for Weierstrass models
    
    y**2 = x**3 + f*x + g
    
    of elliptic curves whose coefficients f and g
    are polynomials of a single variable

    Parameters
    ----------
    fcoeff : coefficients of the polynomial f
    gcoeff : coefficients of the polynomial g
    
    
    Methods
    -------
    get_f : returns the polynomial f
    get_g : returns the polynomial g
    get_D : returns the discriminant D = 4*f**3 + 27*g**2
    get_discriminant_locus : returns the roots of D=0. Roots are sorted
    according to the real part of the locus
    get_number_discriminant_points : returns the number of points in the
    discriminant locus
        
    """
            
    def __init__(self,fcoeff,gcoeff):
        self.u_rot_phase = None
        self.d_min_real = 0.0
        self.d_min = 1.0
        self.f = None
        self.g = None
        self.D = None
        self.disc_locus = None
        u_rot_counter = 0

        while self.d_min_real \
                        < MINIMUM_DISCRIMINANT_LOCUS_SEPARATION * self.d_min \
                and u_rot_counter < U_PLANE_ROTATION_STEPS:
            
            rot_angle = cmath.pi * u_rot_counter / U_PLANE_ROTATION_STEPS
            self.u_rot_phase = np.exp(1j * rot_angle)
            new_fcoeff = rotate_poly(fcoeff, self.u_rot_phase)
            new_gcoeff = rotate_poly(gcoeff, self.u_rot_phase)

            if u_rot_counter > 0:
                logging.info(('\n**************************'
                            +'\nBranch cuts are too close.'
                            +' Will rotate the u-plane by an angle {}\n'
                            +'**************************\n').format(rot_angle))
            
            self.initialize_w_model(new_fcoeff, new_gcoeff)
            u_rot_counter += 1

        self.num = len(self.disc_locus)
        self.lowest_y_coord = min(self.disc_locus.imag)        


    def initialize_w_model(self, fcoeff, gcoeff):
        self.f = np.poly1d(fcoeff)
        self.g = np.poly1d(gcoeff)
        self.D = 4*self.f**3+27*self.g**2
        #Accounting for cancellations of higher order terms in
        #discriminant
        for i, coeff in enumerate(self.D.c):
            if np.absolute(coeff) > NEGLIGIBLE_BOUND:
                self.D = np.poly1d(self.D.c[i:])
                break
        
        temp_disc_locus = np.array(sorted(self.D.r,cmp=real_part))
        self.disc_locus = [temp_disc_locus[0]]
        for entry in temp_disc_locus:
            if np.abs(entry-self.disc_locus[-1]) > NEGLIGIBLE_BOUND:
                self.disc_locus.append(entry)
        self.disc_locus = np.array(self.disc_locus)
        
        distance_half_matrix = [[self.disc_locus[i] - self.disc_locus[j] \
                                for i,x in enumerate(self.disc_locus[:j])] 
                                        for j,y in enumerate(self.disc_locus)]
        distances = [item for sublist in distance_half_matrix \
                                                        for item in sublist]

        self.d_min = min(map(abs, distances))
        self.d_min_real = min(map(abs, map(get_real_part, distances)))
    
                        
    def get_f(self):
        return self.f
        
    def get_g(self):
        return self.g
    
    def get_D(self):
        return self.D
        
    def get_discriminant_locus(self):
        return self.disc_locus 
        
    def get_number_discriminant_points(self):
        return self.num

    


    


class WeierstrassModelWithPaths(WeierstrassModel):
    """
    WeierstrassModelWithPaths(fcoeff,gcoeff,path_data)
    
    is a class for Weierstrass models
    
    y**2 = x**3 + f*x + g
    
    of elliptic curves whose coefficients f and g
    are polynomials of a single variable.

    Parameters
    ----------
    fcoeff : coefficients of the polynomial f
    gcoeff : coefficients of the polynomial g
    path_data : path data for evolving roots given by
    a 2-entry list consisting of a starting point and
    a direction vector, both complex numbers

    Methods
    -------
    get_path_data : returns the path data, path_data
    get_paths : returns an array, paths
    Each entry of paths is a path associated to a
    discriminant locus d of the Weierstrass model.
    A path is a list of seven points whose
    entries are complex numbers denoting a point
    of the segment. Along path[0]~path[2]
    the path is initialized, while along the points
    path[2]~path[7] the path encircles the locus d.
    get_init_rts : returns the array init_rts
    init_rts has three entries, the three roots at the
    initial point of all the paths sorted in ascending
    order of their real part.
    
    Appendix
    --------
    The roots of the RHS of the Weierstrass model are
    sorted in ascenting order of the real part.
    
    The elliptic curve corresponding to the Weierstrass
    model at the initial point of the path
    is such that the branch cuts are chosen to
    cut between the roots 0~1
    (cut 1), and root 2~infinity (cut 2). The
    A-cycle is chosen to encircle cut 1 while, the
    B-cycle is chosen to connect cut 1 and 2. 
    """
            
    def __init__(self,fcoeff,gcoeff,path_data=None):

        WeierstrassModel.__init__(self,fcoeff,gcoeff)
        if path_data == None:
            self.path_data=sample_path_data(self.disc_locus)
        else:
            self.path_data=path_data
        init_p = self.path_data[1]
        
        ### Will declare once and for all the basepoint that is chosen 
        ### to compute monodromies and hence trivialize the lattice 
        ### fibration. All charges are expressed in terms of periods
        ### of the lattice at this point!
        self.base_point = init_p
        
        init_poly = np.poly1d([1,0,self.f(init_p),self.g(init_p)])
        self.init_rts = sorted(init_poly.r,cmp=real_part)        
        self.paths=[]
        self.paths_for_periods = []
        for loc in self.disc_locus:
            self.paths.append(construct_path(loc,self.path_data))

        self.x_rotation_consistency_check = True
        

    
    def get_path_data(self):
        return self.path_data

    def get_paths(self):
        return self.paths        
                        
    def get_init_rts(self):
        return self.init_rts      
  
    def compute_initial_periods(self):
        """
        Computing periods and their derivatives at the 
        base point (called init_p in the __init__ method above).

        Let the three roots e_1, e_2, e_3 be ordered with
        Re(e_1) < Re(e_2) < Re(e_3), then we call gamma_1
        the cycle stretching from e_2 to e_3, and 
        gamma_2 the cycle stretching from e_1 to e_2.
        We have <gamma_1, gamma_2> = 1 in this way.

        Using some carefully engineered functions, keeping track of 
        branch cuts in the x-plane, we compute the period of dx/y 
        along gama_1 and call that eta_0, similarly we denote the
        period along gamma_2 by beta_0.
        We also compute their derivatives.
        """

        u_0 = self.path_data[1]
        init_poly_0 = np.poly1d([1,0,self.f(u_0),self.g(u_0)])
        rts_0 = sorted(init_poly_0.r,cmp=real_part)
        e1_0, e2_0, e3_0 = [complex(rts_0[0]),\
                            complex(rts_0[1]),\
                            complex(rts_0[2])
                            ]

        ### NOTE: eta is related to period B, while beta 
        ### is related to period A. That's because
        ### with current conventions <B, A> = +1
        eta_0 = period_A(e1_0, e2_0, e3_0)[0]
        beta_0 = period_B(e1_0, e2_0, e3_0)[0]

        return [eta_0, beta_0]
      
        
class WeierstrassProto(WeierstrassModelWithPaths):
    """
    A class for Weierstrass models to observe
    monodoromies conveniently.
    """
    def __init__(self,fcoeff,gcoeff,path_data=None):
        WeierstrassModelWithPaths.__init__(self,fcoeff,gcoeff,\
                                            path_data)
        
    def show_path(self,dnum):
        fig=plt.figure()
        rpart, ipart = [], []
        for loc in self.disc_locus:
            rpart.append(loc.real)
            ipart.append(loc.imag)
        plt.plot(rpart,ipart,'ro')
        rpart, ipart = [], []
        for loc in self.paths[dnum]:
            rpart.append(loc.real)
            ipart.append(loc.imag)
        plt.plot(rpart,ipart,'k-')
        plt.show()        
    
    def plot_discriminant(self):
        fig=plt.figure()
        rpart, ipart = [], []
        for loc in self.disc_locus:
            rpart.append(loc.real)
            ipart.append(loc.imag)
        plt.plot(rpart,ipart,'ro')
        plt.show()
        
       
                                 
#### General purpose helper functions ####
                        
def cx_to_coords(cx_number):
    return [cx_number.real,cx_number.imag]
  
      
def angles_from_roots(roots):
    """
    angles_from_roots(roots)
    
    extracts the angles between three
    roots, "roots."
    
    Parameters
    ----------
    
    """
    angles = []
    for i in range(3):
        angles.append(np.angle(roots[(i+1)%3]-roots[i]))
    return angles



#### Functions for constructing paths and contous ####
        
def sample_path_data(dLoc):
    """
    sample_paths(dLoc) computes data for sample paths
    along which the evolution of A and B-cycles are tracked.
    It returns an initial point for the paths and the minimum
    spacing between adjacent loci.
    
    
    Parameters
    ----------
    dLoc: an array of aligned discriminant loci
    
    
    Output
    ------
    returns [a,p,l]
    a : minimum spacing
    p : initial point for path to start from
    l : minimum distance between two roots
    """
    min_y = min([loc.imag for loc in dLoc])
    spacing = 0.45 * min([dLoc[i+1].real-dLoc[i].real \
                                            for i in range(len(dLoc)-1)])
    initPoint = 1j*(min_y-5*spacing)+(dLoc[0].real-2*spacing)
    
    ### A manual choice
    # initPoint = 4.0-1.0j

    ### Another possible automatized choice of the basepoint
    ###
    # max_spacing = max([dLoc[i+1].real-dLoc[i].real \
    #                            for i in range(len(dLoc)-1)])
    # min_spacing = min([dLoc[i+1].real-dLoc[i].real \
    #                            for i in range(len(dLoc)-1)])
    # initPoint = 1j*(min_y - 1.1 * max_spacing) \
    #                            + (dLoc[0].real-0.2 * min_spacing)

    logging.debug('The basepoint for the Wmodel is : {}'.format(initPoint))
    min_len=np.absolute(dLoc[0]-dLoc[1])
    for c in combinations(dLoc,2):
        min_len=min(min_len,np.absolute(c[0]-c[1]))
    return [spacing,initPoint,min_len]            



def construct_path(d,path_data):
    """
    construct_path(d,path_data)
    
    constructs the path for determining the cycle
    type of discriminant locus "d" of a Weierstrass
    model based on the path data.
    
    Parameters
    ----------
    d: discriminant locus (complex number)
    path_data: a list of path data. path_data[0]
    is the "spacing." path_data[1] is the starting
    point of the path. path_data[2] is the minimum
    distance between roots.
    
    Returns
    -------
    path
    
    path is a list of seven points whose
    entries are complex numbers denoting a point
    of the segment. Along path[0]~path[2]
    the path is initialized, while along the points
    path[2]~path[7] the path encircles the locus d.
    """
    
    delta, start, edge = path_data
    edge *= 0.8*np.sqrt(0.5)
    path = []
    
    #Building path to lower left corner of square surrounding
    #discriminant locus
    path.append(start)
    path.append(1j*start.imag+(d-delta).real)
    path.append(d-delta-(1j)*edge)  
    

    #Building path of square surrounding discriminant locus
    #going in counter-clockwise direction
    path.append(d+edge+(-1j)*edge)
    path.append(d+edge+(1j)*edge)
    path.append(d-edge+(1j)*edge)
    path.append(d-edge-(1j)*edge)
    path.append(d-delta-(1j)*edge)  

    return path                


#### Functions related to monodromy ####

def monodromy(G):
    """
    monodromy(G)
    
    computes the monodromy matrix MM
    from the braiding array G
    
    Parameter
    ---------
    G: an array of elements [G_1,G_2,...] of the braid group
    
    Returns
    -------
    returns MM
    MM is a 2x2 matrix representation of the group element
    G_1 G_2 ...
    i.e.,
    MM = rho(G_1).rho(G_2). ....
    """
    def letter_to_matrix(g):
        if g == 'X':
            return np.matrix([[1,0],[1,1]])
        elif g == 'x':
            return np.matrix([[1,0],[-1,1]])
        elif g == 'Y':
            return np.matrix([[1,-1],[0,1]])
        elif g == 'y':
            return np.matrix([[1,1],[0,1]])
        else:
            raise ValueError("SL(2,Z) element "+g+" unrecognized.")
    
    MM=np.matrix([[1,0],[0,1]])
    for g in G:
        MM = np.dot(MM,letter_to_matrix(g))
        
    return MM
        
    
    
def invert_monodromy(MM):
    """
    monodromy_inverse(MM)
    
    computes the inverse of a
    monodromy matrix MM

    Parameter
    ---------
    MM: a monodromy matrix
    
    Returns
    -------
    returns MM^(-1)
    """
    
    if np.shape(MM) != (2,2):
        raise ValueError("Monodromy matrix has wrong dimensions.")
                
    return np.matrix([[MM[1,1],-MM[0,1]],[-MM[1,0],MM[0,0]]])
        
        
      
def monodromy_at_point_via_path(init_root, d, f, g, path_data, w_model,\
                                 ts_ipath=1000, ts_mon=1000, option=None):
    """
    monodromy_at_point_via_path(init_root,d,f,g,path_data,ts_ipath,ts_mon)
    
    yields the monodrompy matrix of locus "d" of a
    Weierstrass model with coefficients f and g,
    based on the open subset of the base manifold
    defined by path_data.
    
    Parameters
    ----------
    init_root: root at initial point
    d: discriminant locus (complex number)
    f, g: Weierstrass coefficients.
    init_data : a list [a,b,c] of three roots
    a, b, c of the Weierstrass model.
    path_data: a list of path data. path_data[0]
    is the "spacing." path_data[1] is the starting
    point of the path. path_data[2] is the minimum
    distance between roots.
    ts_ipath: timesteps for initial path
    ts_mon: timesteps for monodromy
    option: when option='print' the function prints
    the monodromy group elements
            
    Returns
    ------
    The monodromy matrix around the point "d" via
    the path constructed using path_data
    """

    total_path = construct_path(d,path_data)
    initial_path = total_path[:3]
    circ_path = total_path[2:]
    MM = np.matrix([[1,0],[0,1]])

    original_init_root = init_root

    control_var = 'fine'
    
    segments=len(initial_path)-1
    for s in range(segments):
        ev_result=evolve_to_get_braiding(init_root,f,g,\
                                        [initial_path[s],initial_path[s+1]],\
                                        ts_ipath)
        if ev_result == 'rotate':
            control_var = 'rotate'
            break
        
        else:
            if option == 'p':
                print("Initialization path "+str(s)+": "+str(ev_result[0]))
            MM = np.dot(MM,monodromy(ev_result[0]))
            init_root = ev_result[1]
        
    init_MM = MM
    
    if control_var == 'fine':   #proceed with circular part of the path
        segments=len(circ_path)-1
        for s in range(segments):
            ev_result=evolve_to_get_braiding(init_root,f,g,\
                                            [circ_path[s],circ_path[s+1]],\
                                            ts_mon)
            if ev_result == 'rotate':
                control_var = 'rotate'
                break

            else:
                if option == 'p':
                    print("Encircling path "+str(s)+":"+str(ev_result[0]))
                MM = np.dot(MM,monodromy(ev_result[0]))
                init_root = ev_result[1]
        
    if control_var == 'fine':
        mon_matrix = np.dot(MM,invert_monodromy(init_MM))
        return mon_matrix

        ### The following script was needed when we
        ### had monodromies with negative eigenvalues
        ### keep for now.
        ###
        # ### Now we make sure to return a monodromy
        # ### matrix whose eigenvalue is +1, not -1
        # ### this would otherwise cause trouble
        # ### with kwalls undergoing the wrong monodromy, 
        # ### then intersecting each other with negative
        # ### pairing.
        # eigen_vals = list(LA.eig(mon_matrix)[0])
        # if eigen_vals[0] >= 0 and eigen_vals[1] >= 0:
        #     return mon_matrix
        # else:
        #     return (-1 * mon_matrix)

    elif control_var == 'rotate':
        ### Here we trigger a rotation of the x-plane
        ### this will cause the evaluation of ALL 
        ### monodromies to start over from scratch.
        w_model.x_rotation_consistency_check = False
        return np.matrix([[0, 0], [0, 0]])

            
def monodromy_at_point_wmodel(n,wmodel,ts_ipath=1000,ts_mon=1000,option=None):
    """
    monodromy_at_point_wmodel(n,wmodel,ts_ipath,ts_mon,option)
    
    yields the monodrompy matrix of the n'th
    discriminant locus of the Weierstrass model
    wmodel.
    
    Parameters
    ----------
    n: index of discriminant locus 0~(D-1), where D is the
    number of discriminant loci
    ts_ipath: timesteps for initial path
    ts_mon: timesteps for monodromy
    option: when option='p' the function prints
    the monodromy group elements
            
    Returns
    ------
    MM
    
    MM: The monodromy matrix around the n'th discriminant locus
    of the Weierstrass model wmodel
    """
    if n >= wmodel.get_number_discriminant_points():
        raise ValueError("There are less than "+str(n+1)+\
                           " discriminat loci "+\
                           "in the Weierstrass model.")
    return monodromy_at_point_via_path(wmodel.get_init_rts(),\
                                        wmodel.get_discriminant_locus()[n],\
                                        wmodel.get_f(),wmodel.get_g(),\
                                        wmodel.get_path_data(),\
                                        wmodel,\
                                        ts_ipath,ts_mon,option)

                                            

#### Evloution Functions ####                                                
                                                                                                                                            
def evolve_roots_and_angles(init_data,f,g,start_end,timesteps=1000):
    """
    evolve_roots_and_angles(init_data,f,g,start_end,timesteps)
    
    evolves the roots and angles between roots of a given
    Weierstrass model from the start point to end point
    by tracking the roots closely.
    
    This one works!
    
    Parameters
    ----------
    init_data : a list [a,b,c,ab,bc,ca] of three roots
    a, b, c of the Weierstrass model and three angles
    arg(b-a), arg(c-b) and arg(a-c) at the starting point.
    f, g : The coefficient polynomials of the
    Weierstrass models
    start_end : list of two points. start_end[0] is the
    start point, while start_end[1] is the end point.
    
    Returns
    ------
    list1
    
    list1: a list of [re(a),im(a),re(b),im(b),re(c),im(c),ab,bc,ca]
    during the evolution.
    """

    #RHS of the differential equation for evolution of
    #the roots and angles
    def crds_to_roots(crds):
        return [crds[0]+1j*crds[1],crds[2]+1j*crds[3],crds[4]+1j*crds[5]]
    
    def roots_to_crds(roots):
        return [roots[0].real,roots[0].imag,roots[1].real,roots[1].imag,\
                 roots[2].real,roots[2].imag]
    
    def is_degenerate(roots):
        for i in range(3):
            if roots[i] == roots[(i+1)%3]:
                return True
        else:
            return False            
                            
    def sort_roots(rts,last_rts):
        sorted_rts=[None]*3
        max_diff_cons = 0
        for i,a in enumerate(last_rts):
            min_diff = np.absolute(a-rts[0])+1
            for b in rts:
                if np.absolute(a-b) < min_diff:
                    sorted_rts[i] = b
                    min_diff = np.absolute(a-b)
            max_diff_cons = max(max_diff_cons,min_diff)
        if is_degenerate(sorted_rts):
            raise ValueError('Roots not being properly distinguished: '+\
                               'timesteps too sparse.')
        else:            
            return [sorted_rts,max_diff_cons]    
#        return [sorted_rts, max_diff_cons]

    def angle_diff(rts,last_rts):
        angles = []
        for i in range(3):
            new_diff = (rts[(i+1)%3] - rts[i])
            old_diff = (last_rts[(i+1)%3] - last_rts[i])
            angles.append(np.angle(new_diff/old_diff))
        return angles

    def minimum_diff_rts(rts):
        diffs = []
        for i in range(3):
            diff = rts[(i+1)%3]-rts[i]
            diffs.append(np.absolute(diff))
        return min(diffs)

    #setting up initial conditions
    roots_and_angles=[]
    start, end = start_end
    vector=end-start
    min_diff_rts = np.absolute(init_data[0]-init_data[1])
    max_diff_cons = 0
    data = []
    for i in range(3):
        data.append(init_data[i].real)
        data.append(init_data[i].imag)
        diff = np.absolute(init_data[(i+1)%3]-init_data[i])
        min_diff_rts = min(min_diff_rts,diff)
    for i in range(3):
        data.append(init_data[i+3])
    roots_and_angles.append(data)
    time = np.linspace(0,1.,timesteps)
    
    #solve differential equation
    for t in time[1:]:
        poly = np.poly1d([1,0,f(start+t*vector),g(start+t*vector)])
        roots = poly.r
        #print 'Time : '+str(t)
        #print 'Roots : '+str(roots)
        last_roots = crds_to_roots(roots_and_angles[-1][:6])
        roots, cons = sort_roots(roots,last_roots)
        if is_degenerate(roots):
            print 'Last roots: '+str(last_roots)
            print 'Current roots:'+str(roots)
            print 'Path data:'+str(start)+' to '+str(end)
            print 'Location'+str(start+t*vector)
            raise ValueError('Roots not being properly distinguished: '+\
                               'timesteps too sparse.')
        angles = angle_diff(roots,last_roots)
        for k in range(3):
            angles[k] += roots_and_angles[-1][k+6]
        roots_and_angles.append(roots_to_crds(roots) + angles)
        min_diff_rts = min(min_diff_rts,minimum_diff_rts(roots))
        max_diff_cons = max(max_diff_cons,cons)
        
        
    #if difference between roots smaller than timesteps, alert
    if min_diff_rts < max_diff_cons:
        print 'Timesteps too large to be reliable'
                    
    return roots_and_angles



def evolve_to_get_braiding(init_data,f,g,start_end,timesteps=1000):
    """
    get_braiding(init_data,f,g,start_end,timesteps)
    
    evolves the roots and finds the braiding of roots
    given Weierstrass model from the start point to end point
    by tracking the roots closely.
    
    Parameters
    ----------
    init_data : a list [a,b,c] of three roots
    a, b, c of the Weierstrass model at the starting point.
    f, g : The coefficient polynomials of the
    Weierstrass models
    start_end : list of two points. start_end[0] is the
    start point, while start_end[1] is the end point.
    
    Returns
    ------
    [G,last_roots]
    
    G : an element of PSL(2,Z), given by a list of generators
    X, x=X^(-1), Y, y=Y^(-1)
    last_roots : The roots of the Weierstrass model at the
    end point of the path 
    """

    def crds_to_roots(crds):
        return [crds[0]+1j*crds[1],crds[2]+1j*crds[3],crds[4]+1j*crds[5]]
    
    def roots_to_crds(roots):
        return [roots[0].real,roots[0].imag,roots[1].real,roots[1].imag,\
                 roots[2].real,roots[2].imag]
    
    def is_degenerate(roots):
        for i in range(3):
            if roots[i] == roots[(i+1)%3]:
                return True
        else:
            return False            
                            
    def sort_roots(rts,last_rts):
        sorted_rts=[None]*3
        max_diff_cons = 0
        for i,a in enumerate(last_rts):
            min_diff = np.absolute(a-rts[0])+1
            for b in rts:
                if np.absolute(a-b) < min_diff:
                    sorted_rts[i] = b
                    min_diff = np.absolute(a-b)
            max_diff_cons = max(max_diff_cons,min_diff)
        if is_degenerate(sorted_rts):
            # raise ValueError('Roots not being properly distinguished: '+\
            #                    'timesteps too sparse.')
            print '\nsort_roots:: Roots not being properly distinguished: '+\
                       'timesteps too sparse. '+ \
                       '\n\nrts:\n%s\n\nlast_rts:\n%s\n\nsorted_rts\n%s'\
                        % (rts,last_rts,sorted_rts)
            return ['rotate', 'rotate']
        else:            
            return [sorted_rts,max_diff_cons]    
#        return [sorted_rts, max_diff_cons]

    def braiding_action(old_rts,new_rts):
        def indices_left_to_right(rts):
            real_rts = [x.real for x in rts]
            right = max(real_rts)
            left = min(real_rts)
            indices = [-1,-1,-1]
            for i, r in enumerate(real_rts):
                if r == left:
                    indices[0] = i
                elif r == right:
                    indices[2] = i
                else:
                    indices[1] = i
            for i in indices:
                if i == -1:
                    raise ValueError('Cannot sort roots for braiding.')
            return indices
                            
        i_old = indices_left_to_right(old_rts)        
        i_new = indices_left_to_right(new_rts)
        
        if i_old == i_new:
            return 'I'
        elif i_old[0]==i_new[1] and i_old[1]==i_new[0]:
            if old_rts[i_old[0]].imag > old_rts[i_old[1]].imag:
                if new_rts[i_old[0]].imag <= new_rts[i_old[1]].imag:
                    print(str(old_rts))
                    print(str(new_rts))
                    raise ValueError('Cannot distinguish x and X.')
                else:
                    return 'X'
            else:
                if new_rts[i_old[0]].imag > new_rts[i_old[1]].imag:
                    print(str(old_rts))
                    print(str(new_rts))
                    raise ValueError('Cannot distinguish x and X.')
                return 'x'
        elif i_old[1]==i_new[2] and i_old[2]==i_new[1]:
            if old_rts[i_old[1]].imag > old_rts[i_old[2]].imag:
                if new_rts[i_old[1]].imag <= new_rts[i_old[2]].imag:
                    print(str(old_rts))
                    print(str(new_rts))
                    raise ValueError('Cannot distinguish y and Y.')
                else:
                    return 'Y'
            else:
                if new_rts[i_old[1]].imag > new_rts[i_old[2]].imag:
                    print(str(old_rts))
                    print(str(new_rts))
                    raise ValueError('Cannot distinguish y and Y.')
                return 'y'
        else:
            logging.debug('\nbraiding_action: Timesteps too small to extract'\
                                                                + ' braiding.')
            return 'rotate'
                
    def minimum_diff_rts(rts):
        diffs = []
        for i in range(3):
            diff = rts[(i+1)%3]-rts[i]
            diffs.append(np.absolute(diff))
        return min(diffs)


    #setting up initial conditions
    start, end = start_end
    vector=end-start
    min_diff_rts = np.absolute(init_data[0]-init_data[1])
    max_diff_cons = 0
    data = []
    for i in range(3):
        data.append(init_data[i].real)
        data.append(init_data[i].imag)
        diff = np.absolute(init_data[(i+1)%3]-init_data[i])
        min_diff_rts = min(min_diff_rts,diff)
    last_roots = init_data
    time = np.linspace(0,1.,timesteps)
    G = []
    
    control_var = 'fine'

    #solve differential equation
    for t in time[1:]:
        poly = np.poly1d([1,0,f(start+t*vector),g(start+t*vector)])
        roots = poly.r
        #print 'Time : '+str(t)
        #print 'Roots : '+str(roots)
        roots, cons = sort_roots(roots,last_roots)
        # if is_degenerate(roots):
        #     print 'Last roots: '+str(last_roots)
        #     print 'Current roots:'+str(roots)
        #     print 'Path data:'+str(start)+' to '+str(end)
        #     print 'Location'+str(start+t*vector)
        #     raise ValueError('Roots not being properly distinguished: '+\
        #                        'timesteps too sparse.')

        if roots == 'rotate':
            point = start+t*vector
            print '\nTrouble sorting roots:\nu = %s\n' % point
            control_var = 'rotate'
            break

        else:
            b = braiding_action(last_roots,roots)
            
            if b != 'I' and b != 'rotate':
                G.append(b)
            elif b == 'rotate':
                control_var = 'rotate'
                break

            min_diff_rts = min(min_diff_rts,minimum_diff_rts(roots))
            max_diff_cons = max(max_diff_cons,cons)
            last_roots = roots
        
    #if difference between roots smaller than timesteps, alert
    if min_diff_rts < max_diff_cons:
        print 'Timesteps too large to be reliable'
                    
    if control_var == 'fine':
        return [G,last_roots]

    elif control_var == 'rotate':
        return 'rotate'




#Test Functions to test evolution functions

def animate_roots_and_angles_path(wmodel,dnum,xspan=5,yspan=2,\
                            timesteps=1000,steps=10,path=None,breaks=None):
    
    def crds_to_roots(crds):
        return [crds[0]+1j*crds[1],crds[2]+1j*crds[3],crds[4]+1j*crds[5]]
                                        
    def boxspan(list,ratio):
        minval=min(list)
        maxval=max(list)
        return (ratio*minval-maxval,ratio*maxval-minval)

    #extract data from wmodel
    init_rts=wmodel.get_init_rts()
    init_angles=angles_from_roots(init_rts)
    f=wmodel.get_f()
    g=wmodel.get_g()
    dLoc=wmodel.get_discriminant_locus()
    if path == None:
        path = wmodel.get_paths()[dnum]
    dLoc_x, dLoc_y = [], []
    for loc in dLoc:
        dLoc_x.append(loc.real)
        dLoc_y.append(loc.imag)
    space = wmodel.get_path_data()[2]
                
    #set figure size    
    fig = plt.figure()
    init_rt_crds=[cx_to_coords(rt) for rt in init_rts]  
    ax_roots = fig.add_subplot(1, 3, 2, aspect='equal')
    ax_angles = [fig.add_subplot(3, 3, 3, aspect='equal'),\
                 fig.add_subplot(3, 3, 6, aspect='equal'),\
                 fig.add_subplot(3, 3, 9, aspect='equal')]
    ax_point = fig.add_subplot(1, 3, 1,)
    ax_roots.set_xlim(*boxspan([c[0] for c in init_rt_crds],xspan))
    ax_roots.set_ylim(*boxspan([c[1] for c in init_rt_crds],yspan))
    for i in range(3):
        ax_angles[i].set_xlim(-1,1)
        ax_angles[i].set_ylim(-1,1)
    ax_point.set_xlim(min(dLoc_x)-2*space,max(dLoc_x)+2*space)
    ax_point.set_ylim(min(dLoc_y)-2*space,max(dLoc_y)+2*space)
                  
    #initializing lines and setting styles
    lines_roots=[]
    lines_angles=[[],[],[]]
    lines_point=[]
    colors='rgb'
    for i in range(3):
        line = ax_roots.plot([],[],colors[i]+'-')[0]
        lines_roots.append(line)
    for i in range(3):
        line = ax_angles[i].plot([],[],'k-')[0]
        lines_angles[i].append(line)
        line = ax_angles[i].plot([],[],'r-',linewidth=2)[0]
        lines_angles[i].append(line)
        line = ax_angles[i].plot([],[],'ro',markeredgecolor='r')[0]
        lines_angles[i].append(line)
    line=ax_point.plot([],[],'bx')[0]
    lines_point.append(line)
    line=ax_point.plot([],[],'k-')[0]
    lines_point.append(line)
    line=ax_point.plot([],[],'r-',linewidth=2)[0]
    lines_point.append(line)
    line=ax_point.plot([],[],'ro',markeredgecolor='r')[0]
    lines_point.append(line)

    #run numerics, extract x, y coordinates   
    root_xcoords=[[],[],[]]
    root_ycoords=[[],[],[]]
    angle_xcoords=[[],[],[]]
    angle_ycoords=[[],[],[]]
    point_xcoords=[]
    point_ycoords=[]

    segments=len(path)-1
    if breaks!=None:
        segments=min(segments,breaks)
    for s in range(segments):
        ev_roots=evolve_roots_and_angles(init_rts+init_angles,
                                         f,g,[path[s],path[s+1]],timesteps)
        for i, t in enumerate(np.linspace(0,1.,timesteps)):
            point_xcoords.append(((1-t)*path[s]+t*path[s+1]).real)
            point_ycoords.append(((1-t)*path[s]+t*path[s+1]).imag)
            for k in range(3):
                root_xcoords[k].append(ev_roots[i][2*k])
                root_ycoords[k].append(ev_roots[i][2*k+1])
                angle_xcoords[k].append(np.cos(ev_roots[i][k+6]))
                angle_ycoords[k].append(np.sin(ev_roots[i][k+6]))
        init_rts=crds_to_roots(ev_roots[-1][:6])
        init_angles=angles_from_roots(init_rts)    
            
    
    #plt.plot(xcoords[0],ycoords[0],'b-',xcoords[3],ycoords[3],'ro')
    
    # initialization function: plot the background of each frame
    def init():
        for line in lines_roots:
            line.set_data([], [])
        for l in lines_angles:
            for line in l:
                line.set_data([],[])
        for line in lines_point:
            line.set_data([],[])    
        return (lines_roots+lines_angles[0]+lines_angles[1]+lines_angles[2]+\
                 lines_point),

    # animation function.  This is called sequentially
    def animate(i):
        for t, line in enumerate(lines_roots):
            line.set_data(root_xcoords[t][:steps*i],root_ycoords[t][:steps*i])
        for t, l in enumerate(lines_angles):
            l[0].set_data(angle_xcoords[t][:steps*i],\
                                            angle_ycoords[t][:steps*i])
            l[1].set_data(angle_xcoords[t][steps*i-20:steps*i],\
                          angle_ycoords[t][steps*i-20:steps*i])
            l[2].set_data(angle_xcoords[t][steps*i-1],\
                                            angle_ycoords[t][steps*i-1])
        lines_point[0].set_data(dLoc_x,dLoc_y)
        lines_point[1].set_data(point_xcoords[:steps*i],\
                                            point_ycoords[:steps*i])
        lines_point[2].set_data(point_xcoords[steps*i-20:steps*i],\
                                point_ycoords[steps*i-20:steps*i])
        lines_point[3].set_data(point_xcoords[steps*i-1],\
                                            point_ycoords[steps*i-1])
        return (lines_roots+lines_angles[0]+lines_angles[1]+lines_angles[2]+\
                lines_point),

    frame_number = len(root_xcoords[0])/steps

    anim = FuncAnimation(fig, animate, init_func=init,
                              frames=frame_number, interval=10, blit=False)
    plt.show()






#Test Functions to test evolution functions


if __name__=="__main__":
    #wmodel=WeierstrassProto([1+1j,-20j,3+1j],[1j,15,-2],[1j,10])
    
    #Nice Model with six loci
    #wmodel=WeierstrassProto([7+1j,-20j+7,3+1j],[-1j+2,15+2j,-2+12j])
    
    #Seiberg-Witten
    # rot=np.exp(np.pi*1j/10.0 * 0.0)
    # xr=np.exp(np.pi*1j/10.0)
    # wmodel=WeierstrassProto(np.array([-1.0/3*rot**2,0,1.0/4.0])/xr**2,
    #                         np.array([-2.0/27*rot**3,0,1.0/12*rot,0])/xr**3)
                            
    #D3 model
    #wmodel=WeierstrassProto([-1j*1.0,0],[1.0])
    
    #SU(2) Nf=1
    #rot=np.exp(np.pi*1j/3.0)
    #wmodel=WeierstrassProto([-3.0*rot**2,0,2.0],[2.0*rot**3,0,-2.0*rot,1.0])
    
    #Arccos model
    #rot=np.exp(np.pi*1j/8.0)
    #wmodel=WeierstrassProto([-1.0],[rot**2*1.0,0,0])

    #N=2*
    rot =1.0
    xr = 1.0
    m = 3.0 + 0.0j
    g2 = -1.0 + 0.0j
    g3 = 1.0 - 1.0j
    wmodel=WeierstrassProto(\
        np.array([ \
                (- 432.0 * g2 )*(rot**2)/1728.0,\
                (- 324.0 * g3 * (m ** 2))*rot/1728.0, \
                ((-9.0 / 4.0) * (g2 ** 2) * (m ** 4))/1728.0\
                ])/xr**2,
        np.array([\
                ( - 432.0 * g3 )/1728.0*(rot**3),\
                ((- 18.0 * (g2 ** 2) * (m ** 2))/1728.0)*(rot**2),\
                ((- (27.0 / 4.0) * g2 * g3 * (m ** 4))/1728.0)*rot,\
                (((1.0 / 32.0) * (g2 ** 3) * (m ** 6) - (27.0 / 16.0) * (g3 ** 2) * (m ** 6))/1728.0)\
                ])/xr**3)

    

    dnum=0
    animate_roots_and_angles_path(wmodel,dnum,4,3,\
                                 timesteps=5000,steps=120,\
                                 path=None,breaks=None)

    mon1 = monodromy_at_point_wmodel(0,wmodel,5000,5000,option='p')

    # mon2 = monodromy_at_point_wmodel(1,wmodel,5000,5000,option='p')

    print('\nMonodromy at locus 1: ')
    print(str(mon1))
    
    # print('\nMonodromy at point 2: ')
    # print(str(mon2))

    # from elliptic_fibration import monodromy_eigencharge as ME
    # print ME(mon1)
    # print ME(mon2)

