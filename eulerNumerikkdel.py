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
        ax = calculate_acceleration(p, x[i], vx, k)
        vx += ax*dt
        x[i+1] = x[i] + vx*dt
        y[i+1] = get_y(p, x[i+1])
        t[i+1] = t[i] + dt
    # eksempelplot:
    # m = 958
    t_192, x_192, y_192 = parse_track('192-autotrack.txt')
    f, (xplot, yplot) = plt.subplots(2, sharex=True)
    xplot.plot(t, x)
    xplot.plot(t_192, x_192)
    xplot.set_title('x(t)')
    yplot.plot(t, y)
    yplot.plot(t_192, y_192)
    yplot.set_title('y(t)')
    plt.show()


def parse_track(filename):
    t_slow = []
    x_slow = []
    y_slow = []
    with open(filename) as f:
        for _ in range(2):
            next(f)
        i = 0
        for line in f:
            parsedLine = line.replace(',', '.').split('\t')
            t_slow.append(float(parsedLine[0]))
            x_slow.append(float(parsedLine[1]))
            y_slow.append(float(parsedLine[2][:-1]))
            i += 1
    return np.array(t_slow), np.array(x_slow), np.array(y_slow)

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


def calculate_acceleration(p, x, vx, k):
    g = 9.81
    a = g*np.sin(alpha(p, x))*5/7 - k*vx
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


simulate(9.57, 0.01, -0.585, 0, 0.05)
# blå er numerisk, oransje er autotracket
# plot_bane("banetrack.txt", -0.7, 0.7, 0.01)
