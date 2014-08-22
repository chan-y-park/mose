import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Actual Animation Function

def grow_and_animate_trajectories(trajectories):
    fig = plt.figure()
    ax = plt.axes(xlim=(-8, 8), ylim=(-8, 8))
    lines=[]
    for t in trajectories:
        line = ax.plot([],[],t.get_color()+'-')[0]
        lines.append(line)

    # initialization function: plot the background of each frame
    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    # animation function.  This is called sequentially
    def animate(i):
        for t, line in enumerate(lines):
            trajectories[t].evolve()
            line.set_data(trajectories[t].get_xcoordinates(),
                          trajectories[t].get_ycoordinates())
        return lines,

    anim = FuncAnimation(fig, animate, init_func=init,
                              frames=2000, interval=10, blit=False)
    plt.show()


def animate_trajectories(trajectories, steps=10):
    fig = plt.figure()
    ax = plt.axes(xlim=(-8, 8), ylim=(-8, 8))
    lines=[]
    xcoords = []
    ycoords = []
    for t in trajectories:
        line = ax.plot([],[],t.get_color()+'-')[0]
        lines.append(line)
        xcoords.append(t.get_xcoordinates())
        ycoords.append(t.get_ycoordinates())

    # initialization function: plot the background of each frame
    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    # animation function.  This is called sequentially
    def animate(i):
        for t, line in enumerate(lines):
            line.set_data(xcoords[t][:steps*i],ycoords[t][:steps*i])
        return lines,

    anim = FuncAnimation(fig, animate, init_func=init,
                              frames=2000, interval=10, blit=False)
    plt.show()