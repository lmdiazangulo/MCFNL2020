import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

def gauss(x,sigma,t=0,c=1,A=1):
    return A*np.exp(-((x-c*t)/(sigma*np.sqrt(2)))**2)

x = np.linspace(-200, 200, num=int(1001))
t = np.linspace(0, 10, num=int(1001))


# plt.plot(x,gauss(x,10,50))
# plt.show()


f = gauss(x,10,t)


# --- Creates animation ---
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 2, 1)
ax1 = plt.axes(xlim=(x[0], x[-1]), ylim=(-1.1, 1.1))
line1,    = ax1.plot([], [], 'o', markersize=1)
timeText1 = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)


def init():
    line1.set_data([], [])
    timeText1.set_text('')
    return line1, timeText1

def animate(i):
    line1.set_data(x, f[i,:])
    timeText1.set_text('Time = %2.1f' % (t[i]*1e9))
    return line1, timeText1

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(t), interval=50, blit=True)

plt.show()
