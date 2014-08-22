import cmath
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



# Sample Locus and Trajectory Classes
# Trajectory class describes object moving in a
# quartic double-well potential.

class Locus(object):

    def __init__(self, name, boundary_conditions):
        self.name=name
        self.xcoordinates = [boundary_conditions[0][0]]
        self.ycoordinates = [boundary_conditions[0][1]]
        self.xvelocity = boundary_conditions[1][0]
        self.yvelocity = boundary_conditions[1][1]
        

    def __str__(self):
        return 'Trajectory name : '+self.name 

    def get_last_coordinate(self):
        return self.coordinates[-1]


            
class Trajectory(Locus):
    """The trajectory class.

    Attributes: coordinates, periods, degeneracy, phase, charge, parents,\
    boundary_condition, count (shared).
    Methods: evolve, terminate (?).
    Arguments of the instance: (initial_charge,degeneracy,phase,parents,boundary_condition)
    """

    #count = 0

    def __init__(self, name, parents, boundary_conditions, color='b'):
        #self.degeneracy = degeneracy
        #self.phase = phase
        self.name = name
        self.parents = parents
        self.xcoordinates = [boundary_conditions[0][0]]
        self.ycoordinates = [boundary_conditions[0][1]]
        self.xvelocity = boundary_conditions[1][0]
        self.yvelocity = boundary_conditions[1][1]
        #self.initial_charge = initial_charge
        #Trajectory.count += 1
        # Make trajectory evolve automatically at this point?
        self.color = color

    def __str__(self):
        return 'Trajectory name : '+self.name 

    def evolve(self):
        def force(z):
            return (-4*(z**3)+68*z)
            
        x, y = self.xcoordinates[-1], self.ycoordinates[-1]        
        new_coordinates=z+self.velocity*0.001
        self.coordinates.append(new_coordinates)
        self.velocity=self.velocity+force(z)*0.001
        print "Evolving trajectory "+self.name+"."+\
   	       "Current position: "+str(new_coordinates)

		
    def get_last_coordinate(self):
        return self.coordinates[-1]
    
    def get_velocity(self):
        return self.velocity
    
    def get_coordinates(self):
        return self.coordinates
        
    def get_name(self):
        return self.name
    
    def get_color(self):
        return self.color

                
# Sample Data
        
a=Locus("point1",(0,0))
path = Trajectory("path1",(a,),(0,0.5+0.8j))



# Checking Trajectory

def test_trajectory(trajectory):
    for i in range(10000):
        print i
        trajectory.evolve()

    plt.plot(trajectory.get_xcoordinates,trajectory.get_ycoordinates)
    plt.show()

#test_trajectory(path)



# Actual Animation Function

def animate_trajectories(trajectories):
    fig = plt.figure()
    ax = plt.axes(xlim=(-8, 8), ylim=(-8, 8))
    lines=[]
    for t in trajectories:
        line = ax.plot([],[],t.get_color()+'-')
        lines.append(line)

    # initialization function: plot the background of each frame
    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    # animation function.  This is called sequentially
    def animate(i):
        for t in trajectories:
            t.evolve()
            datax=[z.real for z in trajectory.get_xcoordinates()]
            datay=[z.imag for z in trajectory.get_ycoordinates()]
            line.set_data(datax, datay)
        return line,

    anim = FuncAnimation(fig, animate, init_func=init,
                              frames=200, interval=10, blit=False)
    plt.show()


animate_trajectory(path)
