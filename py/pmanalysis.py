#!/usr/bin/python
# pmanalysis.py
# Main driver for the phenotype microarray analysis
#
# Author: Daniel A Cuevas
# Created on 22 Nov. 2013
# Updated on 10 Apr. 2015

from __future__ import absolute_import, division, print_function
import argparse
import PMUtil as util
import PMData
import GrowthCurve
import PMFigures
import operator
import pylab as py
import warnings


###############################################################################
# Utility methods
###############################################################################

def curveFilter(clone, rep, w, curve, pmData):
    '''Determine if growth curve passes filters'''
    ODmax = 0.18  # Maximum optical density in first two hours
    # Only check from 30 minute to 2 hour mark
    if [x for x in curve[1:5] if x >= ODmax]:
        # Set filter to True
        pmData.setFilter(clone, rep, w, True)


def printFiltered(pmData):
    '''Print out filtered data'''
    # Get list of filters -- list of tuples
    # Tuple = (clone, source, condition, replicate, [OD values])
    data = pmData.getFiltered()
    # filter_curve file: curves of filtered samples
    fhFilter = open('{}/filter_curves_{}.txt'.format(outDir, outSuffix), 'w')
    fhFilter.write('sample\treplicate\tmainsource\tsubstrate\twell\t')
    fhFilter.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
    fhFilter.write('\n')
    for tup in data:
        # Unpack tuple and obtain data
        clone, source, w, od = tup
        (ms, gc) = pmData.wells[w]

        # Print sample information
        fhFilter.write('{}\t{}\t{}\t{}\t{}\t'.format(clone, rep, ms, gc, w))
        # Print OD readings
        fhFilter.write('\t'.join(['{:.3f}'.format(x) for x in od]))
        fhFilter.write('\n')
    fhFilter.close()


def curveFit(pmData, c, w, finalDataReps, finalDataMean):
    # Iterate through replicates and determine logistic parameters for each
    first = True
    for rep in pmData.replicates[c]:
        if debugOut and plateFlag:
            (ms, gc) = pmData.wells[w]
            info = '\t'.join(c, rep, ms, gc, w)
            util.printStatus('DEBUG: Processing {}'.format(info))
        elif debugOut:
            info = '\t'.join(c, rep, w)
            util.printStatus('DEBUG: Processing {}.'.format(info))

        # Create GrowthCurve object for sample
        currCurve = pmData.getODCurve(c, w, rep)
        gCurve = GrowthCurve.GrowthCurve(currCurve, pmData.time)

        # Save sample GrowthCurve object for records
        finalDataReps[c][w][rep] = gCurve

        # Add to temporary multi-dim array for calculating mean
        if first:
            tempRepData = py.array([gCurve.y0,
                                    gCurve.asymptote,
                                    gCurve.maxGrowthRate,
                                    gCurve.lag],
                                   ndmin=2)
            first = False
        else:
            add = [gCurve.y0, gCurve.asymptote,
                   gCurve.maxGrowthRate, gCurve.lag]
            tempRepData = py.concatenate((tempRepData, [add]))

        if debugOut:
            msg = 'a={}, mgr={}, lag={}'.format(gCurve.asymptote,
                                                gCurve.maxGrowthRate,
                                                gCurve.lag)
            util.printStatus('DEBUG: parameters for '
                             '{} {} {}: {}'.format(c, rep, w, msg))

    # Create logistic curve from replicates' mean parameters
    meanParams = py.mean(tempRepData, axis=0)
    finalDataMean[c][w]['curve'] = GrowthCurve.logistic(pmData.time,
                                                        *meanParams[0:4])
    # Determine growth level
    gl = GrowthCurve.calcGrowth(finalDataMean[c][w]['curve'],
                                meanParams[1])
    meanParams = py.concatenate((meanParams, [gl]))
    finalDataMean[c][w]['params'] = meanParams
    return (finalDataReps, finalDataMean)

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
parser.add_argument('-f', '--filter', action='store_true',
                    help='Apply filtering to growth curves')
parser.add_argument('-g', '--newgrowth', action='store_true',
                    help='Apply new growth level calculation')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Increase output for status messages')
parser.add_argument('-p', '--plate', action='store_true',
                    help='Input wells are based on a plate')
parser.add_argument('-i', '--images', action='store_true',
                    help='Generate images and graphs')

args = parser.parse_args()
inputFile = args.infile
outSuffix = args.outsuffix if args.outsuffix else 'out'
outDir = args.outdir
filterFlag = args.filter
newGrowthFlag = args.newgrowth
verbose = args.verbose
debugOut = args.debug
plateFlag = args.plate
imagesFlag = args.images

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
    if plateFlag:
        util.printStatus('Found {} samples and {} growth conditions.'.format(
            pmData.numClones, pmData.numConditions))
    else:
        util.printStatus('Plate option not given.')
        util.printStatus('Found {} samples.'.format(pmData.numClones))
if debugOut:
    # Print out number of replicates for each clone
    for c, reps in pmData.replicates.items():
        util.printStatus('DEBUG: {} has {} replicates.'.format(c, len(reps)))


# Perform filter
if filterFlag:
    util.printStatus('Performing filtering...')
    # Iterate through clones
    for c, repDict in pmData.dataHash.items():
        # Iterate through replicates
        for rep, wellDict in repDict.items():
            # Iterate through wells
            for w, odDict in wellDict.items():
                # Perform filter check
                curveFilter(c, rep, w, odDict['od'], pmData)

    util.printStatus('Filtering complete.')
    if verbose:
        util.printStatus('Filtered {} growth '
                         'curves.'.format(pmData.numFiltered))

elif verbose:
    util.printStatus('Filtering option not given -- no filtering performed.')


# Create growth curves and logistic models
util.printStatus('Processing growth curves and creating logistic models...')
finalDataReps = {}
finalDataMean = {}
# Iterate through clones
for c in pmData.clones:
    finalDataReps[c] = {w: {} for w in pmData.wells}
    finalDataMean[c] = {w: {} for w in pmData.wells}

    # Iterate through media sources
    for w in pmData.wells:
        (finalDataReps, finalDataMean) = curveFit(pmData, c, w,
                                                  finalDataReps,
                                                  finalDataMean)

util.printStatus('Processing complete.')


###############################################################################
# Output Files
###############################################################################

util.printStatus('Printing output files...')
# Print out filtered data if set
if filterFlag:
    printFiltered(pmData)

# Print out plate info accordingly
if plateFlag:
    plateInfo = 'mainsource\tsubstrate\tplate\twell'
else:
    plateInfo = 'well'

# logistic_params_sample file: logistic curve parameters for each sample
fhLPSample = open('{}/logistic_params_sample_'
                  '{}.txt'.format(outDir, outSuffix), 'w')
fhLPSample.write('sample\treplicate\t{}\ty0\tlag\t'.format(plateInfo))
fhLPSample.write('maximumgrowthrate\tasymptote\tgrowthlevel\tsse\tgrowth\n')

# logistic_curves_sample file: logistic curves for each sample
fhLCSample = open('{}/logistic_curves_sample_'
                  '{}.txt'.format(outDir, outSuffix), 'w')
fhLCSample.write('sample\treplicate\t{}\t'.format(plateInfo))
fhLCSample.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCSample.write('\n')

# logistic_params_mean file. logistic curve parameters (mean)
fhLPMean = open('{}/logistic_params_mean_'
                '{}.txt'.format(outDir, outSuffix), 'w')
fhLPMean.write('sample\t{}\ty0\tlag\t'.format(plateInfo))
fhLPMean.write('maximumgrowthrate\tasymptote\tgrowthlevel\tgrowth\n')


# logistic_curves_mean file. logistic curve (mean)
fhLCMean = open('{}/logistic_curves_mean_'
                '{}.txt'.format(outDir, outSuffix), 'w')
fhLCMean.write('sample\t{}\t'.format(plateInfo))
fhLCMean.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCMean.write('\n')

# Sort well numbers for print out
if plateFlag:
    ws = [(x[0], int(x[1:])) for x in pmData.wells.keys()]
else:
    ws = [(x[0], int(x[1:])) for x in pmData.wells]
sortW = sorted(ws, key=operator.itemgetter(0, 1))
sortW = ['{}{}'.format(x[0], x[1]) for x in sortW]
# Iterate through clones
for c, wellDict in finalDataReps.items():
    # Iterate through wells
    for w in sortW:
        if plateFlag:
            (ms, gc) = pmData.wells[w]
            pn = pmData.plateName[c]
            pInfo = '\t'.join((ms, gc, pn, w))
        else:
            pInfo = w

        # Process mean information
        try:
            curr = finalDataMean[c][w]
        except KeyError:
            # KeyError occurs when all replicates were filtered out so does not
            # exist in final hash
            continue

        # Print sample information
        fhLPMean.write('{}\t{}\t'.format(c, pInfo))
        fhLCMean.write('{}\t{}\t'.format(c, pInfo))

        # Print parameters
        curve = curr['curve']
        y0 = curr['params'][0]
        asymptote = curr['params'][1]
        mgr = curr['params'][2]
        lag = curr['params'][3]
        gLevel = curr['params'][4]
        growth = GrowthCurve.growthClass(gLevel)
        fhLPMean.write('\t'.join(['{:.3f}'.format(x)
                                  for x in (y0, lag, mgr, asymptote, gLevel)]))
        fhLPMean.write('\t{}'.format(growth))
        fhLPMean.write('\n')

        # Print logistic curve
        fhLCMean.write('\t'.join(['{:.3f}'.format(x)
                                  for x in curve]))
        fhLCMean.write('\n')

        # Iterate through replicates
        for rep in pmData.replicates[c]:
            try:
                curve = wellDict[w][rep]
            except KeyError:
                # KeyError occurs when replicate was filtered out so does not
                # exist in final hash
                continue

            # Print sample information
            fhLPSample.write('{}\t{}\t{}\t'.format(c, rep, pInfo))
            fhLCSample.write('{}\t{}\t{}\t'.format(c, rep, pInfo))

            # Print parameters
            y0 = curve.y0
            lag = curve.lag
            mgr = curve.maxGrowthRate
            asymptote = curve.asymptote
            gLevel = curve.growthLevel
            sse = curve.sse
            growth = GrowthCurve.growthClass(gLevel)
            fhLPSample.write('\t'.join(['{:.3f}'.format(x)
                                        for x in (y0, lag, mgr, asymptote,
                                                  gLevel, sse)]))
            fhLPSample.write('\t{}'.format(growth))
            fhLPSample.write('\n')

            # Print logistic curve
            fhLCSample.write('\t'.join(['{:.3f}'.format(x)
                                        for x in curve.dataLogistic]))
            fhLCSample.write('\n')

# Close output file handles
fhLPSample.close()
fhLCSample.close()
fhLPMean.close()
fhLCMean.close()
util.printStatus('Output files complete')

# Check if images need to be generated
if imagesFlag:
    util.printStatus('Generating growth level heatmap...')
    PMFigures.heatMap(finalDataMean, pmData, sortW, outDir, plateFlag)

    ###########################################################################
    # Each growth curve plot uses as input a dictionary in the format:
    #     { 'sample_name' :
    #         { 'well_1' : [ ordered OD values ],
    #           'well_2' : [ ordered OD values ],
    #           ...
    #         },
    #         ...
    #     }
    # Note: lists of OD values are numpy arrays
    ###########################################################################

    util.printStatus('Generating samples growth curves...')
    plotData = {}  # Use as input into curve plotter
    for c in finalDataReps:  # c = sample
        for w in sortW:  # w = well
            for r in finalDataReps[c][w]:  # r = rep
                curve = finalDataReps[c][w][r].rawcurve
                if r in plotData:
                    plotData[r][w] = curve
                else:
                    plotData[r] = {w: curve}
        outpath = '{}/{}_growthcurves.png'.format(outDir, c)
        title = '{} Growth Curves'.format(c)  # Super title of plot
        PMFigures.curvePlotter(plotData, sortW, pmData.time,
                               outpath, title, None)

    # Plot data for median and average data
    util.printStatus('Generating median and average growth curves...')
    plotData = {'med': {}, 'avg': {}, 'err': {}}
    for c in finalDataReps:
        # Plot data for individual samples
        # Use to plot median and average curves with error bars
        plotIndivData = {'med': {c: {w: {} for w in sortW}},
                         'avg': {c: {w: {} for w in sortW}},
                         'err': {c: {w: {} for w in sortW}}}
        for w in sortW:
            for rIdx, r in enumerate(finalDataReps[c][w]):
                curve = finalDataReps[c][w][r].rawcurve
                if rIdx == 0:
                    # Initialize numpy array for first
                    allarray = py.array([curve], ndmin=2)
                else:
                    allarray = py.concatenate((allarray, [curve]))
            # After getting all replicates,
            # calculate median and standard error
            medcurve = py.median(allarray, axis=0)
            avgcurve = py.mean(allarray, axis=0)
            stderr = py.std(allarray, axis=0) / py.sqrt(len(allarray))
            if c in plotData['med']:
                plotData['med'][c][w] = medcurve
                plotData['avg'][c][w] = avgcurve
            else:
                plotData['med'][c] = {w: medcurve}
                plotData['avg'][c] = {w: avgcurve}
                plotData['avg'][c] = {w: avgcurve}
            plotIndivData['med'][c][w] = medcurve
            plotIndivData['avg'][c][w] = avgcurve
            plotIndivData['err'][c][w] = stderr
        outpath = '{}/{}_median_growthcurves.png'.format(outDir, c)
        title = '{} Growth Curves Using Median'.format(c)
        PMFigures.curvePlotter(plotIndivData['med'], sortW, pmData.time,
                               outpath, title, plotIndivData['err'])
        outpath = '{}/{}_average_growthcurves.png'.format(outDir, c)
        title = '{} Growth Curves Using Average'.format(c)
        PMFigures.curvePlotter(plotIndivData['avg'], sortW, pmData.time,
                               outpath, title, plotIndivData['err'])
    # Only plot all median and average curves when there is
    # more than 1 sample
    if pmData.numClones > 1:
        outpath = '{}/all_median_growthcurves'.format(outDir)
        title = 'Median Growth Curves'
        PMFigures.curvePlotter(plotData['med'], sortW, pmData.time,
                               outpath, title, None)
        outpath = '{}/all_average_growthcurves'.format(outDir)
        title = 'Average Growth Curves'
        PMFigures.curvePlotter(plotData['avg'], sortW, pmData.time,
                               outpath, title, None)

util.printStatus('Printing complete.')
util.printStatus('Analysis complete.')
