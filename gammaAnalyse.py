import numpy as np
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
    estimates = np.zeros(len(data_table))
    variances = np.zeros(len(data_table))
    all_g = np.zeros(sum([len(data) - 1 for data in data_table]))
    all_index = 0
    for i in range(len(data_table)):
        wip = np.zeros(len(data_table[i])-1)
        for j in range(len(data_table[i])-1):
            pre = data_table[i][j]
            post = data_table[i][j+1]
            wip[j] = (pre.y - post.y)/(pre.y * (post.t - pre.t))
            all_g[all_index] = wip[j]
            all_index += 1
        estimates[i] = np.average(wip)
        variances[i] = np.var(wip, ddof=1)
    return estimates, variances, all_g


def refined_estimate(estimates, variances):
    refined_average = np.average(estimates)
    refined_variance = np.sum(variances)/(variances.size**2)
    return refined_average, refined_variance


def sd(var):
    return np.sqrt(var)


def stderr(var, n):
    return np.sqrt(var/n)


def plotData(dataset):
    plt.plot(dataset, 'ro')
    plt.show()

data = create_data_table('data', ['P1120186.txt', 'P1120191.txt', 'P1120192.txt', '183.txt', '184.txt', '185.txt',
                                  'ONE.txt', 'TWO.txt', 'THREE.txt', 'FOUR.txt'])
estimates, variances, all_gammas = gamma_estimates(data)
print('Averages:', estimates, '\nVariances:', variances)
refined_avg, refined_var = refined_estimate(estimates, variances)
print('Refined avg:', refined_avg, '\nRefined var:', refined_var, '\nStandardavvik: ', sd(refined_var))
standardfeil = stderr(refined_var, len(data))
print('Standard error:', standardfeil)

plotData(all_gammas)
