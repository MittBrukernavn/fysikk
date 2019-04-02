import matplotlib.pyplot as plt
import numpy as np
from math import ceil


def simulate(time_interval, dt, x0, k=0.0):
    filename = "banetrack.txt"
    p = iptrack(filename)
    n = ceil(time_interval/dt)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    t = np.zeros(n + 1)
    vx = np.zeros(n + 1)
    vy = np.zeros(n + 1)
    f = np.zeros(n + 1)
    N = np.zeros(n + 1)
    x[0] = x0
    y[0] = get_y(p, x[0])
    vx[0] = 0
    f[0], N[0] = f_and_n(x[0], vx[0], base_a(p, x[0]), p)
    for i in range(n):
        ax = calculate_acceleration(p, x[i], vx[i], k)
        vx[i+1] = vx[i] + ax*dt
        x[i+1] = x[i] + vx[i + 1]*dt
        y[i+1] = get_y(p, x[i+1])
        t[i+1] = t[i] + dt
        vy[i + 1] = compute_vy(p, vx[i + 1], x[i + 1])
        f[i+1], N[i+1] = f_and_n(x[i+1], vx[i+1], base_a(p, x[i+1]), p)
    # eksempelplot:
    # m = 958
    t_192, x_192, y_192, vx_192, vy_192 = parse_track('192-autotrack.txt')
    labelsize = 15
    titlesize = 20
    func_i_dont_use, (xplot, yplot) = plt.subplots(2)
    xplot.plot(t, x)
    xplot.plot(t_192, x_192)
    xplot.legend(['Numerisk', 'Eksperimentell'], loc='upper right')
    xplot.set_ylabel('x / m', fontsize=labelsize)
    yplot.set_xlabel('t / s', fontsize=labelsize)
    xplot.set_title('x(t)', fontsize=titlesize)
    yplot.plot(t, y)
    yplot.plot(t_192, y_192)
    yplot.legend(['Numerisk', 'Eksperimentell'], loc='upper right')
    yplot.set_xlabel('t / s', fontsize=labelsize)
    yplot.set_ylabel('y / m', fontsize=labelsize)
    yplot.set_title('y(t)', fontsize=titlesize)
    plt.show()
    gunc_i_dont_use, (vxplot, vyplot) = plt.subplots(2)
    vxplot.set_title('vx(t)', fontsize=titlesize)
    vxplot.set_xlabel('t / s', fontsize=labelsize)
    vxplot.set_ylabel('vx / (m/s)', fontsize=labelsize)
    vxplot.plot(t, vx)
    vxplot.plot(t_192, vx_192)
    vxplot.legend(['Numerisk', 'Eksperimentell'], loc='upper right')
    vyplot.set_title('vy(t)', fontsize=titlesize)
    vyplot.set_xlabel('t / s', fontsize=labelsize)
    vyplot.set_ylabel('vy / (m/s)', fontsize=labelsize)
    vyplot.plot(t, vy)
    vyplot.plot(t_192, vy_192)
    vyplot.legend(['Numerisk', 'Eksperimentell'], loc='upper right')
    plt.show()
    plt.plot(x, f)
    plt.title('Friksjonsakselerasjon', fontsize=titlesize)
    plt.xlabel('x-posisjon / m', fontsize=labelsize)
    plt.ylabel('(f/m) / (m/s^2) ', fontsize=labelsize)
    plt.show()
    plt.plot(x, N)
    plt.title('Normalakselerasjon', fontsize=titlesize)
    plt.xlabel('x-posisjon / m', fontsize=labelsize)
    plt.ylabel('(N/m) / (m/s^2)', fontsize=labelsize)
    plt.show()


def f_and_n(x, v, a, p):
    return 9.81*np.sin(alpha(p, x)) - a, ((v**2)/R(p, x)) + 9.81*np.cos(x)


def compute_vy(p, vx, x):
    dp = np.polyder(p)
    dydx = np.polyval(dp, x)
    return vx*dydx


def parse_track(filename):
    t_slow = []
    x_slow = []
    y_slow = []
    vx_slow = []
    vy_slow = []
    with open(filename) as f:
        for _ in range(2):
            next(f)
        for line in f:
            parsedLine = line.replace(',', '.').split('\t')
            t, x, y = float(parsedLine[0]), float(parsedLine[1]), float(parsedLine[2][:-1])
            t_slow.append(t)
            x_slow.append(x)
            y_slow.append(y)
            if len(t_slow) > 1:
                vx_slow.append((x_slow[-1] - x_slow[-2])/(t_slow[-1] - t_slow[-2]))
                vy_slow.append((y_slow[-1] - y_slow[-2]) / (t_slow[-1] - t_slow[-2]))
    vx_slow.append((x_slow[-1] - x_slow[-2]) / (t_slow[-1] - t_slow[-2]))
    vy_slow.append((y_slow[-1] - y_slow[-2]) / (t_slow[-1] - t_slow[-2]))
    return np.array(t_slow), np.array(x_slow), np.array(y_slow), np.array(vx_slow), np.array(vy_slow)


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


def R(p, x):
    return trvalues(p, x)[4]


def iptrack(filename):
    data = np.loadtxt(filename,skiprows=2)
    return np.polyfit(data[:,1],data[:,2],15)


def calculate_acceleration(p, x, vx, k):
    a = base_a(p, x) - k*vx
    return a


def base_a(p, x):
    return 9.81*np.sin(alpha(p, x)) * (5/7)


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
    plt.ylabel('y / m')
    plt.xlabel('x / m')
    plt.show()


simulate(9.57, 0.01, -0.585, 0.05)
# blå er numerisk, oransje er autotracket
# plot_bane("banetrack.txt", -0.7, 0.7, 0.01)
