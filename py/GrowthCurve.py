# GrowthCurve.py
# Bacteria growth curve class to extract parameters
# Use with Phenotype MicroArray analysis pipeline
#
# Author: Daniel A Cuevas
# Created on 21 Nov. 2013
# Updated on 15 Jan. 2014


import pylab as py
import Models


class GrowthCurve:
    '''Bacteria growth curve class'''
    def __init__(self, data, time):
        # data format: multidimensional numpy array
        #              Each inner array is an array of OD values
        #              ordered by time.
        #              This is important for determining the median

        self.dataReps = data  # OD data values (replicates implied)
        self.dataMed = py.median(self.dataReps, axis=0)
        self.time = time  # time values
        self.asymptote = self.__calcAsymptote()
        self.maxGrowthRate, self.mgrTime = self.__calcMGR()
        self.dataLogistic, self.lag = self.__calcLag()
        self.growthLevel = self.__calcGrowth()

    def __calcAsymptote(self):
        '''Obtain the value of the highest OD reading'''
        # Calculate asymptote using a sliding window of 3 data points
        stop = len(self.time) - 3
        maxA = -1
        for idx in xrange(1, stop):
            av = py.mean(self.dataMed[idx:idx + 3])
            if av > maxA:
                maxA = av
        return maxA

    def __calcMGR(self):
        '''Obtain the value of the max growth'''
        # Calculate max growth rate using a sliding window of 4 data points
        stop = len(self.time) - 4
        maxGR = 0
        for idx in xrange(1, stop):

            # Growth rate calculation:
            # (log(i+3) - log(i)) / (time(i+3) - time(i))
            gr = ((py.log(self.dataMed[idx + 3]) - py.log(self.dataMed[idx])) /
                  (self.time[idx + 3] - self.time[idx]))
            if idx == 1 or gr > maxGR:
                maxGR = gr
                t = self.time[idx + 2]  # Midpoint time value

        return maxGR, t

    def __calcLag(self):
        '''Obtain the value of the lag phase using best fit model'''
        logisticData, lag, sseF = Models.Models(self.dataMed, self.dataMed[1],
                                                self.maxGrowthRate,
                                                self.asymptote,
                                                self.time).Logistic()
        return logisticData, lag

    def __calcGrowth(self):
        '''Calculate growth level using an adjusted harmonic mean'''
        return len(self.dataLogistic) / py.sum([(1 / (x + self.asymptote))
                                                for x in self.dataLogistic])
