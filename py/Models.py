# Models.py
# Bacteria growth curve models
#
# Author: Daniel A Cuevas
# Created on 21 Nov. 2013
# Updated on 15 Jan. 2014

import pylab as py


class Models:
    '''Class containing growth curve models using given growth parameters'''
    def __init__(self, data, startOD, maxgrowth, asymptote, time):
        self.data = data
        self.startOD = startOD
        self.maxgrowth = maxgrowth
        self.asymptote = asymptote
        self.time = time

    def Logistic(self):
        '''Create logistic model from data'''
        tStep = self.time[1] - self.time[0]

        # Time vector for calculating lag phase
        timevec = py.arange(self.time[0], self.time[-1], tStep / 2)

        # Try using to find logistic model with optimal lag phase
        # y = p2 + (A-p2) / (1 + exp(( (um/A) * (L-t) ) + 2))
        sse = 0
        sseF = 0

        # Attempt to use every possible value in the time vector as the lag
        # Choose lag that creates best-fit model
        for idx, lag in enumerate(timevec):
            logDataTemp = [self.startOD + ((self.asymptote - self.startOD) /
                                           (1 + py.exp((((self.maxgrowth /
                                                          self.asymptote) *
                                                         (lag - t)) + 2)
                                                       ))
                                           ) for t in self.time]
            sse = py.sum([((self.data[i] - logDataTemp[i]) ** 2)
                          for i in xrange(len(self.data) - 1)])
            if idx == 0 or sse < sseF:
                logisticData = logDataTemp
                lagF = lag
                sseF = sse

        return logisticData, lagF, sseF


# Gompertz model not available/not used right now

#    def Gompertz(self):
#        '''Create Gompertz model from data'''
#        tStep = self.time[1] - self.time[0]
#        # Only go up to time of inflection point
#        timevec = arange(time[0], mGrowTime, tStep)
#
#        # Try using to find Gompertz model with optimal lag phase
#        # y = (A-p2) * exp(-exp( ((um*e/A-p2) * (L-t)) + 1 ))
#        for idx, lag in enumerate(timevec):
#            gompDataTemp = [((self.asymptote - self.startOD) *
#                             exp(-exp((self.maxgrowth * e /
#                                       (self.asymptote - self.startOD) *
#                                      (lag - t)) + 1)))
#                           for t in self.time]
#            sse = sum([(data[i] - gompDataTemp[i]) ** 2
#                       for i in xrange(len(data) - 1])
#            if idx == 0 or sse < sseF:
#                gompertzData = gompDataTemp
#                lagF = lag
#                sseF = sse
#
#        return (gompertData, lagF, sseF)
