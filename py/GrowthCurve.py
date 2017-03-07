# GrowthCurve.py
# Bacteria growth curve class to extract parameters
# Use with Phenotype MicroArray analysis pipeline
#
# Author: Daniel A Cuevas
# Created on 21 Nov 2013
# Updated on 07 Mar 2017


from __future__ import absolute_import, division, print_function
import PMUtil as util
import pylab as py
import scipy.optimize as optimize
from numpy import trapz


class GrowthCurve:
    """Bacteria growth curve class"""
    def __init__(self, data, newGrowth=0):
        # data format
        #   Pandas GroupBy group object
        #   Indices: sample rep well time
        #   Value: OD reading

        self.rawcurve = data.values  # OD values
        self.time = py.array(data.index.get_level_values("time"))
        # Lowest y0 is chosen
        # Instances of condensation issue at beginning of experiment can cause
        # high OD values
        self.y0 = py.amin(self.rawcurve[0:3])
        self.asymptote = self.__calcAsymptote()
        self.maxGrowthRate, self.mgrTime = self.__calcMGR()
        self.y0, self.asymptote, self.maxGrowthRate, self.lag = (
            self.__calcParameters(
                (self.y0, self.asymptote, self.maxGrowthRate, 0.01),
                self.time, self.rawcurve)
        )

        self.dataLogistic = logistic(self.time,
                                     self.y0,
                                     self.asymptote,
                                     self.maxGrowthRate,
                                     self.lag)
        if newGrowth == 2:
            self.growthLevel = calcNewGrowth2(self.dataLogistic,
                                              self.asymptote,
                                              self.y0)
        elif newGrowth == 1:
            self.growthLevel = calcNewGrowth(self.dataLogistic,
                                             self.asymptote,
                                             self.y0)
        else:
            self.growthLevel = calcGrowth(self.dataLogistic, self.asymptote)
        self.glScaled = calcGrowth2(self.dataLogistic, self.asymptote)
        self.expGrowth = calcExpGrowth(self.maxGrowthRate, self.asymptote)

        self.auc_raw = calcAUCData(self.rawcurve, self.time)
        self.auc_rshift = calcShiftAUC(self.auc_raw, self.y0, self.time[-1])
        self.auc_log = calcAUC(self.rawcurve, self.y0, self.lag,
                               self.maxGrowthRate, self.asymptote, self.time)
        self.auc_lshift = calcShiftAUC(self.auc_log, self.y0, self.time[-1])

        self.growthClass = growthClass(self.growthLevel)
        self.sse = sum((self.dataLogistic - self.rawcurve) ** 2)
        self.mse = self.sse / len(self.time)

    def __calcParameters(self, y0, t, raw):
        """Perform curve-fitting optimization to obtain parameters"""
        try:
            results = optimize.minimize(self.__logisticSSE, y0, args=(t, raw),
                                        bounds=((0, None),
                                                (0.01, y0[1] + 0.5),
                                                (0, None),
                                                (0, None)))
        except RuntimeError as e:
            util.printStatus(e)
            util.printStatus(self.rawcurve)
            util.exitScript()

        if not results.success:
            util.printStatus("*" * 55)
            util.printStatus("CurveFit Unsuccessful")
            util.printStatus(results.message)
            util.printStatus("*" * 55)

        return results.x

    def __logisticSSE(self, params, t, y):
        y0, a, mgr, l = params
        return py.sum((logistic(t, y0, a, mgr, l) - y) ** 2)

    def __calcAsymptote(self):
        """Obtain the initial guess of the highest OD reading"""
        # Calculate asymptote using a sliding window of 3 data points
        stop = len(self.time) - 3
        maxA = -1
        for idx in range(1, stop):
            av = py.mean(self.rawcurve[idx:idx + 3])
            if av > maxA:
                maxA = av
        return maxA

    def __calcMGR(self):
        """Obtain the initial guess of the max growth"""
        # Calculate max growth rate using a sliding window of 4 data points
        stop = len(self.time) - 3
        maxGR = 0
        t = 0
        for idx in range(0, stop):
            # Growth rate calculation:
            # (log(i+3) - log(i)) / (time(i+3) - time(i))

            # Be sure to check for zeros
            ya = self.rawcurve[idx]
            yb = self.rawcurve[idx + 2]
            if ya <= 0:
                ya = 0.001
            if yb <= 0:
                yb = 0.001
            gr = ((py.log(yb) - py.log(ya)) /
                  (self.time[idx + 2] - self.time[idx]))
            if idx == 1 or gr > maxGR:
                maxGR = gr
                t = self.time[idx + 1]  # Midpoint time value

        return maxGR, t


def calcGrowth(logistic, asym):
    """
    Calculate growth level using an adjusted harmonic mean
    using a logistic model and its asymptote
    """
    return len(logistic) / py.sum((1 / (logistic + asym)))


def calcNewGrowth(logistic, asym, y0):
    """
    Calculate growth level using an adjusted harmonic mean
    using a logistic model, its asymptote, and its starting OD value
    """
    diff = asym - y0
    return len(logistic) / py.sum((1 / (logistic + diff)))


def calcNewGrowth2(logistic, asym, y0):
    """
    Calculate growth level using an adjusted harmonic mean
    using a logistic model, its asymptote, and its starting OD value
    """
    amplitude = asym - y0
    diff = logistic - y0
    return len(logistic) / py.sum((1 / (amplitude + diff)))


def calcGrowth2(logistic, asym):
    """
    Calculate growth level using an adjusted harmonic mean
    using a logistic model and its asymptote
    """
    return (len(logistic) / (asym * py.sum((1 / (logistic + asym))))) - 1


def calcExpGrowth(mgr, asym):
    """Calculate exponential growth value 'r'"""
    return 4 * mgr / asym


def calcAUC(data, y0, lag, mgr, asym, time):
    """
    Calculate the area under the curve of the logistic function
    using its integrated formula
    [ A( [A-y0] log[ exp( [4m(l-t)/A]+2 )+1 ]) / 4m ] + At
    """

    # First check that max growth rate is not zero
    # If so, calculate using the data instead of the equation
    if mgr == 0:
        auc = calcAUCData(data, time)
    else:
        timeS = time[0]
        timeE = time[-1]
        t1 = asym - y0
        t2_s = py.log(py.exp((4 * mgr * (lag - timeS) / asym) + 2) + 1)
        t2_e = py.log(py.exp((4 * mgr * (lag - timeE) / asym) + 2) + 1)
        t3 = 4 * mgr
        t4_s = asym * timeS
        t4_e = asym * timeE

        start = (asym * (t1 * t2_s) / t3) + t4_s
        end = (asym * (t1 * t2_e) / t3) + t4_e
        auc = end - start

    if py.absolute(auc) == float('Inf'):
        x = py.diff(time)
        auc = py.sum(x * data[1:])
    return auc


def calcAUCData(data, time):
    """Calculate the area under the curve based on data values"""
    ld = len(data)
    lt = len(time)
    if ld != lt:
        raise ValueError("len(data) {} != len(time) {}".format(ld, lt))
    try:
        return trapz(data, time)
    except Exception:
        raise


def calcShiftAUC(auc, y0, tF):
    """
    Calculate the area under the curve minus
    the area below the starting OD
    """
    return auc - (y0 * tF)


def growthClass(gLevel):
    """Determine growth class based on growth level"""
    if gLevel >= 0.75:
        return "+++"
    elif gLevel >= 0.50:
        return "++"
    elif gLevel >= 0.35:
        return "+"
    elif gLevel >= 0.25:
        return "-"
    else:
        return "--"


def logistic(t, y0, a, mgr, l):
    """Logistic modeling"""
    startOD = y0
    exponent = ((4 * mgr / a) * (l - t)) + 2
    try:
        denom = 1 + py.exp(exponent)
        lg = startOD + ((a - startOD) / denom)
    except RuntimeWarning as rw:
        util.printStatus("*" * 55)
        util.printStatus("RuntimeWarning: {}".format(rw))
        util.printStatus("   Exponent value for logistic "
                         "equation is {:.3f}".format(exponent[0]))
        util.printStatus("   This produces a large value in the denominator "
                         "of the logistic equation")
        util.printStatus("   Now setting logistic value to y0: "
                         "{:.3f}".format(startOD))
        util.printStatus("*" * 55)
        lg = startOD

    return lg
