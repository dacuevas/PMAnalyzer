# GrowthCurve.py
# Bacteria growth curve class to extract parameters
# Use with Phenotype MicroArray analysis pipeline
#
# Author: Daniel A Cuevas
# Created on 21 Nov. 2013
# Updated on 16 Jan. 2015


from __future__ import absolute_import, division, print_function
import pylab as py
import scipy.optimize as optimize
import sys


class GrowthCurve:
    '''Bacteria growth curve class'''
    def __init__(self, data, time):
        # data format:  numpy array if OD values

        self.rawcurve = data  # OD data values (replicates implied)
        self.time = time  # time values
        self.y0 = self.rawcurve[1]
        self.asymptote = self.__calcAsymptote()
        self.maxGrowthRate, self.mgrTime = self.__calcMGR()
        self.asymptote, self.maxGrowthRate, self.lag = self.__calcParameters(
            (self.asymptote, self.maxGrowthRate, 0.5), self.time, self.rawcurve)

        self.dataLogistic = logistic(self.time, self.y0,
                                     self.asymptote, self.maxGrowthRate, self.lag)
        self.growthLevel = calcGrowth(self.dataLogistic, self.asymptote)
        self.sse = sum((self.dataLogistic - self.rawcurve) ** 2)

    def __calcParameters(self, y0, t, raw):
        '''Perform curve-fitting optimization to obtain parameters'''
        try:
            results = optimize.minimize(self.__logisticSSE, y0, args=(t, raw),
                                        bounds=((0.01, y0[0]),
                                                (0, None),
                                                (0, None)))
        except RuntimeError as e:
            print(e)
            print(self.rawcurve)
            sys.exit(1)

        return results.x

    def __logisticSSE(self, params, t, y):
        a, mgr, l = params
        return py.sum((logistic(t, self.y0, a, mgr, l) - y) ** 2)

    def __calcAsymptote(self):
        '''Obtain the initial guess of the highest OD reading'''
        # Calculate asymptote using a sliding window of 3 data points
        stop = len(self.time) - 3
        maxA = -1
        for idx in range(1, stop):
            av = py.mean(self.rawcurve[idx:idx + 3])
            if av > maxA:
                maxA = av
        return maxA

    def __calcMGR(self):
        '''Obtain the initial guess of the max growth'''
        # Calculate max growth rate using a sliding window of 4 data points
        stop = len(self.time) - 4
        maxGR = 0
        t = 0
        for idx in range(1, stop):

            # Growth rate calculation:
            # (log(i+3) - log(i)) / (time(i+3) - time(i))
            gr = ((py.log(self.rawcurve[idx + 3]) - py.log(self.rawcurve[idx])) /
                  (self.time[idx + 3] - self.time[idx]))
            if idx == 1 or gr > maxGR:
                maxGR = gr
                t = self.time[idx + 2]  # Midpoint time value

        return maxGR, t


def calcGrowth(logistic, asym):
    '''Calculate growth level using an adjusted harmonic mean
    using a logistic model and its asymptote'''
    return len(logistic) / py.sum((1 / (logistic + asym)))


def growthClass(gLevel):
    '''Determine growth class based on growth level'''
    if gLevel >= 0.75:
        return '+++'
    elif gLevel >= 0.50:
        return '++'
    elif gLevel >= 0.35:
        return '+'
    elif gLevel >= 0.25:
        return '-'
    else:
        return '--'


def logistic(t, y0, a, mgr, l):
    '''Logistic modeling'''
    startOD = y0
    exponent = ((mgr / a) * (l - t)) + 2
    lg = startOD + ((a - startOD) /
                    (1 + py.exp(exponent)))
    return lg
