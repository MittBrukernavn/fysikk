import matplotlib.pyplot as plt
import numpy as np
from math import ceil


def simulate(time_interval, dt, x0, v0):
    filename = "track.txt"
    p = iptrack(filename)
    n = ceil(time_interval/dt)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    t = np.zeros(n + 1)
    vx = v0
    x[0] = x0
    for i in range(n):
        ax = calculate_acceleration(p, x[i])
        vx += ax*dt
        x[i+1] = x[i] + vx*dt
        y[i+1] = get_y(p, x[i])
        t[i+1] = t[i] + dt
    plt.plot(t, x)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.grid()
    plt.show()


def trvalues(p, x):
    y = np.polyval(p, x)
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    ddp = np.polyder(dp)
    d2ydx2 = np.polyval(ddp, x)
    alpha = np.arctan(-dydx)
    R = (1.0+dydx**2)**1.5/d2ydx2
    return [y, dydx, d2ydx2, alpha, R]


def get_y(p, x):
    return np.polyval(p, x)


def alpha(p, x):
    return trvalues(p, x)[3]


def iptrack(filename):
    data = np.loadtxt(filename,skiprows=2)
    return np.polyfit(data[:,1],data[:,2],15)


def calculate_acceleration(p, x):
    g = 9.81
    a = g*np.sin(alpha(p, x))*5/7
    ax = a*np.sin(alpha(p, x))
    return ax

simulate(100,0.01,1,0)