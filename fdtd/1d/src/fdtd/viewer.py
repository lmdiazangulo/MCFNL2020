import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class Animator:

    def __init__(self, mesh, probe,analytical=None, layer=None,dispAndFree=False,fps=5):
            
        ids = probe["indices"]
        gridE = mesh.pos[ids[0]:ids[1]]

        probeTime = probe["time"][:]
        values    = probe["values"][:]



        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(1, 1, 1)
        #upperlim = np.min([np.nanmax(values[1:]),1])
        #lowerlim = np.max([np.nanmin(values[1:]),-1])
        ax1 = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1, 1))
        ax1.grid(color='gray', linestyle='--', linewidth=.2)
        ax1.set_xlabel('X coordinate [m]')
        ax1.set_ylabel('Field')
        line1,    = ax1.plot([], [], '-', markersize=1)
        
        if dispAndFree:
            valuesFree    = probe["valuesFree"][:]
        line2,    = ax1.plot([], [], '-', markersize=1)     #For dispersive and free space
 

        line3,    = ax1.plot([], [], '-', markersize=1)        #For analytical sol
        timeText1 = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

        if layer is not None:
            ax1.fill_betweenx([-1,1], layer.coords[0], layer.coords[-1],color='orange',alpha=0.50)
        def init():
            line1.set_data([ ], [ ])
            if dispAndFree:
                line2.set_data([ ], [ ])
            if analytical:
                line3.set_data(gridE, analytical[i][:])
            timeText1.set_text('')
            return line1,line2, timeText1

        def animate(i):
            line1.set_data(gridE, values[i][:])
            if dispAndFree:
                line2.set_data(gridE, valuesFree[i][:])
            if analytical:
                line3.set_data(gridE, analytical[i][:])
            timeText1.set_text('Time = %.1E [ns]' % (probeTime[i]*1e9))
            return line1,line2, timeText1
            
        animation.FuncAnimation(fig, animate, init_func=init,
            frames=len(probeTime), interval=1000/fps, blit=True)

        plt.show()


