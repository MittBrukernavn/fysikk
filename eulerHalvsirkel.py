import matplotlib.pyplot as plt
import numpy
from math import ceil


def run_eulers_method(time_interval, dt, x0, v0):
    n = ceil(time_interval/dt)
    x = numpy.zeros(n+1)
    t = numpy.zeros(n + 1)
    v = v0
    x[0] = x0
    for i in range(n):
        a = calculate_acceleration(x[i], v)
        v += a*dt
        x[i+1] = x[i] + v*dt
        t[i+1] = t[i] + dt
    plt.plot(t, x)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.grid()
    plt.show()


def calculate_acceleration(x, v):
    k = 0.05
    lk = 0.05
    return - k * x - lk * v

run_eulers_method(100,0.01,1,0)