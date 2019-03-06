import matplotlib.pyplot as plt
import numpy as np
from math import ceil


def simulate(time_interval, dt, x0, v0, k=0.0):
    filename = "banetrack.txt"
    p = iptrack(filename)
    n = ceil(time_interval/dt)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    t = np.zeros(n + 1)
    vx = v0
    x[0] = x0
    y[0] = get_y(p, x[0])
    for i in range(n):
        ax = calculate_acceleration(p, x[i])
        vx += ax*dt - k*vx*dt
        x[i+1] = x[i] + vx*dt
        y[i+1] = get_y(p, x[i+1])
        t[i+1] = t[i] + dt
    plt.plot(t, x)
    # eksempelplot:
    m = 958
    t_192 = np.zeros(m)
    x_192 = np.zeros(m)
    y_192 = np.zeros(m)
    with open('192-autotrack.txt') as f:
        for _ in range(2):
            next(f)
        i = 0
        for line in f:
            parsedLine = line.replace(',', '.').split('\t')
            t_192[i], x_192[i], y_192[i] = float(parsedLine[0]), float(parsedLine[1]), float(parsedLine[2][:-1])
            i += 1
    plt.plot(t_192, x_192)
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
    return a


def plot_bane(filename, minx, maxx, step):
    n = ceil((maxx-minx)/step)
    p = iptrack(filename)
    y = np.zeros(n+1)
    x = np.zeros(n+1)
    x[0] = minx
    y[0] = np.polyval(p, x[0])
    for i in range(n):
        x[i+1] = x[i] + step
        y[i+1] = get_y(p, x[i+1])
    plt.plot(x, y)
    plt.show()


simulate(10, 0.01, -0.585, 0, 0.1)
# blå er numerisk, oransje er autotracket
# plot_bane("banetrack.txt", -0.7, 0.7, 0.01)
