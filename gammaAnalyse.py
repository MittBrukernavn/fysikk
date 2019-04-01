import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

class DataPoint:
    def __init__(self, y, t):
        self.y = y
        self.t = t


def create_data_table(folder, fileNames):
    data_table = []
    for filename in fileNames:
        points = []
        with open(folder+'/'+filename) as f:
            for line in f.readlines()[2:]:
                t, x, y = line.replace(',', '.').split('\t')
                points.append(DataPoint(float(y), float(t)))
        data_table.append(points)
    return data_table


def gamma_estimates(data_table):
    n = sum([len(data) - 1 for data in data_table])
    all_g = np.zeros(n)
    all_t = np.zeros(n)
    all_index = 0
    for i in range(len(data_table)):
        for j in range(len(data_table[i])-1):
            pre = data_table[i][j]
            post = data_table[i][j+1]
            all_g[all_index] = (pre.y - post.y)/(pre.y * (post.t - pre.t))
            all_t[all_index] = pre.t
            all_index += 1
    return all_g, all_t


def refined_estimate(estimates):
    refined_average = np.average(estimates)
    refined_variance = np.var(estimates, ddof=1)
    return refined_average, refined_variance


def stderr(var, n):
    return np.sqrt(var/n)


def plotData(dataset):
    plt.plot(dataset, 'ro')
    plt.ylim(0, 0.2)
    plt.xlabel('Måling nr')
    plt.ylabel('Målt verdi / (1/s)')
    plt.show()


def plot_data_over_time(gammas, times):
    slope, intercept, r_value, p_value, std_err = stats.linregress(times, gammas)
    print("Stigning:", slope, "\nSkjæringspunkt:", intercept, "\nStderr:", std_err,
          "\nP-value given null hypothesis that slope is zero:", p_value, "\nCorr:", r_value)
    regtimes = np.array([0, 10])
    regvals = np.array([intercept, intercept+10*slope])
    plt.plot(times, gammas, 'ro')
    plt.plot(regtimes, regvals)
    plt.ylim(0, 0.2)
    plt.xlabel('Tid / s')
    plt.ylabel('Målt verdi / (1/s)')
    plt.show()


data = create_data_table('data', ['P1120186.txt', 'P1120191.txt', 'P1120192.txt', '183.txt', '184.txt', '185.txt',
                                  'ONE.txt', 'TWO.txt', 'THREE.txt', 'FOUR.txt'])
all_gammas, all_times = gamma_estimates(data)

refined_avg, refined_var = refined_estimate(all_gammas)
print('Refined avg:', refined_avg, '\nRefined var:', refined_var, '\nStandardavvik: ', np.sqrt(refined_var))
standardfeil = stderr(refined_var, len(all_gammas))
print('Standard error:', standardfeil)

plotData(all_gammas)
plot_data_over_time(all_gammas, all_times)
