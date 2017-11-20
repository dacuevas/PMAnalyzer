#!/usr/bin/python
# pmanalysis.py
# Main driver for the phenotype microarray analysis
#
# Author: Daniel A Cuevas
# Created on 22 Nov 2013
# Updated on 14 Nov 2017

from __future__ import absolute_import, division, print_function
import argparse
import PMUtil as util
import PMData
import GrowthCurve
import pylab as py
import warnings
import pandas as pd


###############################################################################
# Utility methods
###############################################################################

def buildLogParams(data):
    """Build DataFrame to hold logistic parameters"""
    d = {"y0": 0.0, "maxgrowth": 0.0, "asymptote": 0.0,
         "lag": 0.0, "growthlevel": 0.0, "glscaled": 0.0,
         "r": 0.0, "auc_raw": 0.0, "auc_rshift": 0.0,
         "auc_log": 0.0, "auc_lshift": 0.0,
         "growthclass": "", "err_mse": 0.0}
    if plateFlag:
        d["mainsource"] = "ms"
        d["compound"] = "cp"

    # Remove time level
    logP_DF = pd.DataFrame(d, index=data.index.droplevel("time"))
    # Reset index to drop duplicated rows caused by time points
    logP_DF.reset_index(inplace=True)
    logP_DF.drop_duplicates(inplace=True)
    logP_DF.set_index(["sample", "rep", "well"], inplace=True)
    return logP_DF


def curveFit(group, args):
    """Perform logistic fit on each growth curve using Pandas GroupBy groups"""
    sample, rep, well = group.name
    dataLogParams = args[0]
    plateFlag = args[1]
    growthFlag = args[2]

    # Perform logistic fitting
    try:
        if plateFlag:
            gCurve = GrowthCurve.GrowthCurve(group["od"], sample, rep, well,
                                             growthFlag)
        else:
            gCurve = GrowthCurve.GrowthCurve(group, sample, rep, well,
                                             growthFlag)
    except Exception as e:
        util.printStatus("sample: {}, rep: {}, well: {}".format(*group.name))
        util.printStatus(e)
        util.exitScript()

    # Add logistic parameters to DataFrame
    time = group.index.get_level_values("time")
    dataLogParams["y0"].loc[sample, rep, well] = gCurve.y0
    dataLogParams["maxgrowth"].loc[sample, rep, well] = gCurve.maxGrowthRate
    dataLogParams["asymptote"].loc[sample, rep, well] = gCurve.asymptote
    dataLogParams["lag"].loc[sample, rep, well] = gCurve.lag
    dataLogParams["growthlevel"].loc[sample, rep, well] = gCurve.growthLevel
    dataLogParams["glscaled"].loc[sample, rep, well] = gCurve.glScaled
    dataLogParams["r"].loc[sample, rep, well] = gCurve.expGrowth
    dataLogParams["auc_raw"].loc[sample, rep, well] = gCurve.auc_raw
    dataLogParams["auc_rshift"].loc[sample, rep, well] = gCurve.auc_rshift
    dataLogParams["auc_log"].loc[sample, rep, well] = gCurve.auc_log
    dataLogParams["auc_lshift"].loc[sample, rep, well] = gCurve.auc_lshift
    dataLogParams["growthclass"].loc[sample, rep, well] = gCurve.growthClass
    dataLogParams["err_mse"].loc[sample, rep, well] = gCurve.mse
    if plateFlag:
        try:
            m = group["mainsource"][0]
        except TypeError as e:
            stars = "*" * 55
            err = "\n" + stars + "\nERROR: trying to obtain 'mainsource' from "
            err += " ".join((sample, rep, well)) + "\n" + stars
            util.printStatus(err)
            util.printStatus(e)
            util.printStatus(group["mainsource"].values)
            util.printStatus(group)
            util.exitScript()
        try:
            c = group["compound"][0]
        except TypeError as e:
            stars = "*" * 55
            err = "\n" + stars + "\nERROR: trying to obtain 'compound' from "
            err += " ".join((sample, rep, well)) + "\n" + stars
            util.printStatus(err)
            util.printStatus(e)
            util.printStatus(group["compound"].values)
            util.printStatus(group)
            util.exitScript()
        try:
            dataLogParams["mainsource"].loc[sample, rep, well] = m
            dataLogParams["compound"].loc[sample, rep, well] = c
        except TypeError as e:
            stars = "*" * 55
            err = "\n" + stars + "\nERROR: trying to use loc() method on "
            err += " ".join((sample, rep, well)) + "\n" + stars
            util.printStatus(err)
            util.printStatus(e)
            util.printStatus(group)
            util.exitScript()

    # Create Series object of logistic values to return
    idxNames = ["sample", "rep", "well", "time"]
    index = zip([group.name] * len(time), time)
    index = [x[0] + tuple([x[1]]) for x in index]
    index = pd.MultiIndex.from_tuples(index, names=idxNames)

    # Create DataFrame object with new logistic curve
    d = {"od": gCurve.dataLogistic}
    if plateFlag:
        d["mainsource"] = m
        d["compound"] = c
    df = pd.DataFrame(d, index=py.round_(time, decimals=3))
    df.index.name = "time"
    return df


def getLogCurve(group, args):
    """Return logistic curve of each Pandas GroupBy group"""
    plateFlag = args[0]
    finalTime = args[1][group.name]  # key = (sample, well)
    time = py.linspace(0.0, finalTime, 100)
    logistic = GrowthCurve.logistic(time,
                                    group["y0"][0],
                                    group["asymptote"][0],
                                    group["maxgrowth"][0],
                                    group["lag"][0])

    # Create DataFrame object with logistic curve
    d = {"od": logistic}
    if plateFlag:
        d["mainsource"] = group["mainsource"][0]
        d["compound"] = group["compound"][0]

    df = pd.DataFrame(d, index=py.round_(time, decimals=3))
    df.index.name = "time"
    return df


###############################################################################
# Argument Parsing
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='Input PM file')
parser.add_argument('outdir',
                    help='Directory to store output files')
parser.add_argument('-o', '--outsuffix',
                    help='Suffix appended to output files. Default is "out"')
parser.add_argument('--debug', action='store_true',
                    help='Output for debugging purposes')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Increase output for status messages')
parser.add_argument('-p', '--plate', action='store_true',
                    help='Input wells are based on a plate')
parser.add_argument('-i', '--images', action='store_true',
                    help='Generate images and graphs')
parser.add_argument('-g', '--growth', type=int,
                    choices=[1, 3], metavar="[1, 3]",
                    help='Use alternative growth level calculations')

args = parser.parse_args()
inputFile = args.infile
outSuffix = args.outsuffix if args.outsuffix else 'out'
outDir = args.outdir
verbose = args.verbose
debugOut = args.debug
plateFlag = args.plate
imagesFlag = args.images
growthFlag = args.growth if args.growth else 0

# Set warnings filter for catching RuntimeWarnings
warnings.filterwarnings('error')

###############################################################################
# Data Processing
###############################################################################

# Parse data file
util.printStatus('Parsing input file...')
pmData = PMData.PMData(inputFile, plateFlag)
util.printStatus('Parsing complete.')
if verbose:
    if not plateFlag:
        util.printStatus('Plate option not given.')
    util.printStatus('Found {} samples and {} wells.'.format(
        pmData.getNumSamples(), pmData.getNumWells()))

# Perform the logistic fitting on the growth curves
# Create a GroupBy Pandas data structure and apply the curvefit
# function to each group
# Grouping is based on sample, replicate, and well. Each group consists of all
# time points and OD values for a single replicate on a single well
# The the growth curve parameters are stored in the dataLogParams DataFrame
util.printStatus('Processing growth curves and creating logistic models...')
dataLogParams = buildLogParams(pmData.DF)
if plateFlag:
    dataLogistic = pmData.DF[["od", "mainsource", "compound"]].groupby(
        level=["sample", "rep", "well"])
else:
    dataLogistic = pmData.DF["od"].groupby(level=["sample", "rep", "well"])
dataLogistic = dataLogistic.apply(curveFit, args=(dataLogParams,
                                                  plateFlag,
                                                  growthFlag))

# Get final time points
finalTimePoints = {}
for name, group in pmData.DF.groupby(level=["sample", "well"]):
    # Add tuple name=(sample, well) to final time points
    t = group.index.get_level_values("time")[-1]
    finalTimePoints[name] = t

# Calculate average parameters for each sample then
# create new logistic curve for each sample
meanParams = dataLogParams.mean(level=["sample", "well"])
if plateFlag:
    # Add mainsource and compounds to group
    #leftMerge = meanParams.reset_index(level=[0], drop=True)
    leftMerge = meanParams.reset_index()
    rightMerge = dataLogParams[["mainsource", "compound"]].reset_index(
        level=[1], drop=True).reset_index().drop_duplicates()
    meanParams = pd.merge(
        leftMerge,
        rightMerge,
        on=["sample", "well"],
        how="left"
    ).set_index(["sample", "well"])

meanParams.drop(pd.Index(["err_mse"]), axis=1, inplace=True)
meanLogCurves = meanParams.groupby(level=["sample", "well"])
meanLogCurves = meanLogCurves.apply(getLogCurve, args=(plateFlag,
                                                       finalTimePoints))
meanParams.loc[:, "growthclass"] = ""

# Calculate mean and median values for printing
dataMean = pmData.getMeanCurves()
dataMedian = pmData.getMedianCurves()

# Calculate new metrics for mean parameters
for name, group in dataLogistic["od"].groupby(level=["sample", "well"]):
    s, w = name
    A = meanParams.loc[s, w]["asymptote"]
    mgr = meanParams.loc[s, w]["maxgrowth"]
    y0 = meanParams.loc[s, w]["y0"]
    lag = meanParams.loc[s, w]["lag"]
    time = group.index.get_level_values("time").unique().values

    # Get mean raw curve to calculate mean auc_raw and auc_rshift
    rc = dataMean.loc[(s, w), "od"].values
    auc_raw = GrowthCurve.calcAUCData(rc, time)
    auc_rshift = GrowthCurve.calcShiftAUC(auc_raw, y0, time[-1])

    # Calculate logistic based on mean values
    log = GrowthCurve.logistic(time, y0, A, mgr, lag)

    if growthFlag == 0:
        gl = GrowthCurve.default_growth(log, A, y0)
    elif growthFlag == 1:
        gl = GrowthCurve.calcNewGrowth(log, A, y0)
    elif growthFlag == 2:
        gl = GrowthCurve.calcGrowth(log, A)
    else:
        gl = GrowthCurve.calcGrowthScore(A, mgr)
    glScaled = GrowthCurve.calcGrowth2(log, A)
    r = GrowthCurve.calcExpGrowth(mgr, A)
    try:
        auc_log = GrowthCurve.calcAUC(log, y0, lag, mgr, A, time)
    except Exception as e:
        util.printStatus("*" * 55)
        util.printStatus("ERROR")
        util.printStatus("*" * 55)
        util.printStatus("sample: {}, well: {}".format(*name))
        util.printStatus(e)
        util.exitScript()
    auc_lshift = GrowthCurve.calcShiftAUC(auc_log, y0, time[-1])
    gClass = GrowthCurve.growthClass(gl)

    cols = ("growthlevel", "glscaled", "r",
            "auc_raw", "auc_rshift", "auc_log", "auc_lshift",
            "growthclass")
    meanParams.loc[(s, w), cols] = (gl, glScaled, r,
                                    auc_raw, auc_rshift,
                                    auc_log, auc_lshift,
                                    gClass)

util.printStatus('Processing complete.')


###############################################################################
# Output Files
###############################################################################

util.printStatus('Printing output files...')

plateCol = []  # Will add the mainsource and compound columns if they exist
if plateFlag:
    plateCol = ["mainsource", "compound"]

# all_median_curves file: median curves for each sample
curvesToPrint = dataMedian.reset_index().set_index(["sample", "well"])
col = plateCol + ["time", "od"]
curvesToPrint.to_csv(
    "{}/all_curves_median_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)


# all_average_curves file: average curves for each sample
curvesToPrint = dataMean.reset_index().set_index(["sample", "well"])
col = plateCol + ["time", "od"]
curvesToPrint.to_csv(
    "{}/all_curves_mean_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)

# logistic_params_sample file: logistic curve parameters for each sample
col = plateCol + ["y0", "lag", "maxgrowth", "asymptote",
                  "growthlevel", "glscaled", "r", "auc_raw", "auc_rshift",
                  "auc_log", "auc_lshift", "growthclass", "err_mse"]
dataLogParams.to_csv(
    "{}/logistic_params_sample_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)

# logistic_curves_sample file: logistic curves for each sample
# Build new logistic curve consisting of 100 data points
# Old method:
#curvesToPrint = dataLogistic.reset_index().set_index(["sample", "rep", "well"])
numPoints = 100
logcurvesdf = pd.DataFrame({"sample": [], "rep": [], "well": [],
                            "mainsource": [], "compound": [],
                            "time": [], "od": []})
for name, group in dataLogParams.groupby(level=["sample", "rep", "well"]):
    finalTime = finalTimePoints[name[0], name[2]]  # key = (sample, well)
    time = py.linspace(0.0, finalTime, numPoints)
    logistic = GrowthCurve.logistic(time,
                                    group["y0"][0],
                                    group["asymptote"][0],
                                    group["maxgrowth"][0],
                                    group["lag"][0])

    tmp_dict = {"sample": name[0],
                "rep": name[1],
                "well": name[2],
                "time": time,
                "od": logistic}
    if plateFlag:
        tmp_dict["mainsource"] = group["mainsource"][0]
        tmp_dict["compound"] = group["compound"][0]
    logcurvesdf = logcurvesdf.append(pd.DataFrame(tmp_dict))

curvesToPrint = logcurvesdf.set_index(["sample", "rep", "well"])
col = plateCol + ["time", "od"]
curvesToPrint.to_csv(
    "{}/logistic_curves_sample_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)

# logistic_params_mean file. logistic curve parameters (mean)
col = plateCol + ["y0", "lag", "maxgrowth", "asymptote",
                  "growthlevel", "glscaled", "r", "auc_raw", "auc_rshift",
                  "auc_log", "auc_lshift", "growthclass"]
meanParams.to_csv("{}/logistic_params_mean_{}.txt".format(outDir, outSuffix),
                  sep="\t", header=True, index=True, float_format="%.3f",
                  columns=col)

# logistic_curves_mean file. logistic curve (mean)
curvesToPrint = meanLogCurves.reset_index().set_index(["sample", "well"])
col = plateCol + ["time", "od"]
curvesToPrint.to_csv(
    "{}/logistic_curves_mean_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)


util.printStatus('Output files complete')
util.printStatus('Analysis complete.')
