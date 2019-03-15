import numpy as np


class DataPoint:
    def __init__(self, y, t):
        self.y = y
        self.t = t


def create_data_table(fileNames):
    data_table = []
    for filename in fileNames:
        points = []
        with open(filename) as f:
            for line in f.readlines()[2:]:
                t, x, y = line.replace(',', '.').split('\t')
                points.append(DataPoint(float(y), float(t)))
        data_table.append(points)
    return data_table


def gamma_estimates(data_table):
    estimates = np.zeros(len(data_table))
    variances = np.zeros(len(data_table))
    for i in range(len(data_table)):
        wip = np.zeros(len(data_table[i])-1)
        for j in range(len(data_table[i])-1):
            pre = data_table[i][j]
            post = data_table[i][j+1]
            wip[j] = (pre.y - post.y)/(pre.y *(post.t - pre.t))
        estimates[i] = np.average(wip)
        variances[i] = np.var(wip, ddof=1)
    return estimates, variances

def refined_estimate(estimates, variances):
    refined_average = np.average(estimates)
    refined_variance = np.sum(variances)/(variances.size**2)
    return refined_average, refined_variance


def sd(var):
    return np.sqrt(var)

def stderr(var, n):
    return np.sqrt(var/n)

data = create_data_table(['P1120186.txt', 'P1120191.txt', 'P1120192.txt'])
estimates, variances = gamma_estimates(data)
print('Averages:', estimates, '\nVariances:', variances)
refined_avg, refined_var = refined_estimate(estimates, variances)
print('Refined avg:', refined_avg, '\nRefined var:', refined_var)
n = 10*len(data)
standardfeil = stderr(refined_var, n)
print('Standard error:', standardfeil)