import cmath
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import odeint
from itertools import combinations



NEGLIGIBLE_BOUND = 0.1**12

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
    so that the 
    get_number_discriminant_points : returns the number of points in the
    discriminant locus
        
    """
            
    def __init__(self,fcoeff,gcoeff):
        def real_part(x,y):
            if x.real>y.real:
                return 1
            else: return -1
            
        self.f = np.poly1d(fcoeff)
        self.g = np.poly1d(gcoeff)
        self.D = 4*self.f**3+27*self.g**2
        #Accounting for cancellations of higher order terms in
        #discriminant
        for i, coeff in enumerate(self.D.c):
            if np.absolute(coeff) > NEGLIGIBLE_BOUND:
                self.D = np.poly1d(self.D.c[i:])
                break
        self.disc_locus = np.array(sorted(self.D.r,cmp=real_part))
        self.num = len(self.disc_locus)
        self.lowest_y_coord = min(self.disc_locus.imag)        
                        
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
    


class WeierstrassModelWithLabels(WeierstrassModel):
    """
    WeierstrassModelWithLabels(fcoeff,gcoeff,cycle_choice,path_data)
    
    is a class for Weierstrass models
    
    y**2 = x**3 + f*x + g
    
    of elliptic curves whose coefficients f and g
    are polynomials of a single variable.

    Parameters
    ----------
    fcoeff : coefficients of the polynomial f
    gcoeff : coefficients of the polynomial g
    cycle_choice : an array given by a permutation
    of 0,1,2 which specifies the choice of an A-cycle
    and a B-cycle. See appendix.
    path_data : path data for evolving roots given by
    a 2-entry list consisting of a starting point and
    a direction vector, both complex numbers
    
    
    Appendix
    --------
    The elliptic curve corresponding to the Weierstrass
    model at the initial point of the path
    is such that the branch cuts are chosen to
    cut between the roots cycle_type[0]~cycle_type[1]
    (cut 1), and cycle_type[2]~infinity (cut 2). The
    A-cycle is chosen to encircle cut 1 while, the
    B-cycle is chosen to connect cut 1 and 2. 
    """
            
    def __init__(self,fcoeff,gcoeff,cycle_choice,path_data=None):
        WeierstrassModel.__init__(self,fcoeff,gcoeff)
        self.cycle_choice=cycle_choice
        if path_data == None:
            self.path_data=sample_path_data(self.disc_locus)
        else:
            self.path_data=path_data
        init_p = self.path_data[1]
        init_poly = np.poly1d([1,0,self.f(init_p),self.g(init_p)])
        self.init_rts = [None]*3
        for i in range(3):
            self.init_rts[i]=init_poly.r[cycle_choice[i]]  
        

    def get_cycle_choice(self):
        return self.cycle_choice

    def get_path_data(self):
        return self.path_data
        
    def get_init_rts(self):
        return self.init_rts
        
        
        
class WeierstrassProto(WeierstrassModelWithLabels):
    """
    A class for Weierstrass models to observe
    monodoromies conveniently.
    """
    def __init__(self,fcoeff,gcoeff,cycle_choice,path_data=None):
        WeierstrassModelWithLabels.__init__(self,fcoeff,gcoeff,\
                                            cycle_choice,path_data)
        self.paths=[]
        for loc in self.disc_locus:
            self.paths.append(construct_path(loc,self.path_data))
    
    def get_paths(self):
        return self.paths
        
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
    spacing = 0.45*min([dLoc[i+1].real-dLoc[i].real for i in range(len(dLoc)-1)])
    initPoint = 1j*(min_y-5*spacing)+(dLoc[0].real-2*spacing)
    min_len=np.absolute(dLoc[0]-dLoc[1])
    for c in combinations(dLoc,2):
        min_len=min(min_len,np.absolute(c[0]-c[1]))
    return [spacing,initPoint,min_len]
    
    
    
def locus_type(d,f,g,init_data,path_data):
    """
    cycle_type(d,f,g,init_data,path_data)
    
    yields the 7-brane type of locus "d" of a
    Weierstrass model with coefficients f and g,
    based on the initial conditions init_data
    on the open subset of the base manifold
    defined by path_data.
    
    Parameters
    ----------
    d: discriminant locus (complex number)
    f, g: Weierstrass coefficients.
    init_data : a list [a,b,c] of three roots
    a, b, c of the Weierstrass model.
    path_data: a list of path data. path_data[0]
    is the "spacing." path_data[1] is the starting
    point of the path. path_data[2] is the minimum
    distance between roots.
            
    Returns
    ------
    The type of 
    """
    pass



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
    [re(a),im(a),re(b),im(b),re(c),im(c),ab,bc,ca]
    
    at the end point of evolution.
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
        min_diff_cons = np.absolute(rts[0]-last_rts[0])+1
        for i,a in enumerate(last_rts):
            min_diff = np.absolute(a-rts[0])+1
            for b in rts:
                if np.absolute(a-b) < min_diff:
                    sorted_rts[i] = b
                    min_diff = np.absolute(a-b)
                    min_diff_cons = min(min_diff_cons,min_diff)
#        if is_degenerate(sorted_rts):
#            raise ValueError('Roots not being properly distinguished: '+\
#                               'timesteps too sparse.')
#        else:            
#            return [sorted_rts,min_diff_cons]    
        return [sorted_rts, min_diff_cons]

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
    min_diff_cons = min_diff_rts
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
        min_diff_cons = min(min_diff_cons,cons)
        
        
    #if difference between roots smaller than timesteps, alert
    if min_diff_rts < min_diff_cons:
        print 'Timesteps too large to be reliable'
                    
    return roots_and_angles







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
            l[0].set_data(angle_xcoords[t][:steps*i],angle_ycoords[t][:steps*i])
            l[1].set_data(angle_xcoords[t][steps*i-20:steps*i],\
                          angle_ycoords[t][steps*i-20:steps*i])
            l[2].set_data(angle_xcoords[t][steps*i-1],angle_ycoords[t][steps*i-1])
        lines_point[0].set_data(dLoc_x,dLoc_y)
        lines_point[1].set_data(point_xcoords[:steps*i],point_ycoords[:steps*i])
        lines_point[2].set_data(point_xcoords[steps*i-20:steps*i],\
                                point_ycoords[steps*i-20:steps*i])
        lines_point[3].set_data(point_xcoords[steps*i-1],point_ycoords[steps*i-1])
        return (lines_roots+lines_angles[0]+lines_angles[1]+lines_angles[2]+\
                lines_point),

    frame_number = len(root_xcoords[0])/steps

    anim = FuncAnimation(fig, animate, init_func=init,
                              frames=frame_number, interval=10, blit=False)
    plt.show()





if __name__=="__main__":
    #wmodel=WeierstrassProto([1+1j,-20j,3+1j],[1j,15,-2],\
    #                        [0,1,2],[1j,10])
    wmodel=WeierstrassProto([7+1j,-20j+7,3+1j],[-1j+2,15+2j,-2+12j],\
                            [0,1,2],None)
    #wmodel=WeierstrassProto([-1.0/3,0,1.0/4.0],[-2.0/27,0,1.0/12,0],[0,1,2],None)
    animate_roots_and_angles_path(wmodel,1,2,3,\
                                  timesteps=5000,steps=120,\
                                  path=None,breaks=None)
