#!/usr/bin/python
# pmanalysis.py
# Main driver for the phenotype microarray analysis
#
# Author: Daniel A Cuevas
# Created on 22 Nov. 2013
# Updated on 10 Nov. 2014

import argparse
import sys
import time
import datetime
import PMData
import GrowthCurve
import operator
import pylab as py
import matplotlib.pyplot as plt
import pprint


###############################################################################
# Utility methods
###############################################################################

def timeStamp():
    '''Return time stamp'''
    t = time.time()
    fmt = '[%Y-%m-%d %H:%M:%S]'
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def printStatus(msg):
    '''Print status message'''
    print('{} {}'.format(timeStamp(), msg), file=sys.stderr)
    sys.stderr.flush()


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


def printHeatMap(data, clones, wells, outDir, plateInfo=False):
    '''Make growth heatmap after curve fitting and analysis'''
    #finalDataMean[c][w]['params'] = meanParams
    #        else:
    #            retArray = py.concatenate((retArray,
    #                                       py.array([currCurve])))
    if plateInfo:
        newWells = [x[1] for x in [pmData.wells['{}{}'.format(w[0], w[1])] for w in wells]]
    else:
        newWells = ['{}{}'.format(w[0], w[1]) for w in wells]
    first = True
    for clone in clones:
        tmpArr = []
        for well in wells:
            w = '{}{}'.format(well[0], well[1])
            tmpArr.append(data[clone][w]['params'][4])

        if first:
            plotData = py.array(tmpArr, ndmin=2)
            first = False
        else:
            plotData = py.concatenate((plotData, [tmpArr]))


    ######################################################
    # Plotting
    ######################################################
    numClones = len(clones)
    numWells = len(newWells)

    # Width is 15 inches
    # All measurements are in inches
    width = 15;

    # Height determined by number of clones
    # Cap at 12 inches
    height = len(clones)
    if height > 12:
        height = 12

    # Fontsize of x axis determined by
    # number of wells on x axsis
    # 13 was determined by trial and error
    if numWells > 13:
        xfontsize = numWells / 13
    else:
        xfontsize = 10
    yfontsize = 10

    # Create figure and axis
    fig, ax = plt.subplots()
    fig.set_size_inches(width,height)

    # Create heatmap object using pcolor
    hm = ax.pcolor(plotData, cmap=plt.cm.Greys, edgecolor='black', vmin=0, vmax=1.5)

    # Create color bar legend
    mini = py.amin(plotData)
    maxi = py.amax(plotData)
    cbticks = [mini, maxi, 0.25, 0.75]
    cblabs = ['min', 'max', 'no growth', 'growth']
    cbticks, cblabs = zip(*sorted(zip(cbticks, cblabs)))
    cbar = fig.colorbar(hm, orientation='horizontal')
    cbar.set_ticks(cbticks)
    cbar.set_ticklabels(cblabs)


    # Move x axis to top
    ax.xaxis.tick_top()
    ax.yaxis.tick_left()

    # Align x and y tick marks to center of cells
    ax.set_xticks(py.arange(0, numWells)+0.5)
    ax.set_yticks(py.arange(0, numClones)+0.5)

    # Set tick labels
    ax.set_xticklabels(labels=newWells, minor=False, rotation=90, fontsize=xfontsize)
    ax.set_yticklabels(labels=clones, minor=False)
    ax.axis('tight')
    # Remove tick marks lines
    plt.tick_params(axis='both', left='off', right='off', bottom='off', top='off')

    plt.savefig('{}/growthlevels.png'.format(outDir), dpi=100, bbox_inches='tight')


def curvePlot(data, wells, time):
    '''Plot growth curve plots on multi-faceted figure. Separate graphs based on well.
    Plot clone replicates together.'''
    f, axarr = plt.subplots(8, 12)
    colors = ['r', 'b', 'g']
    for idx, w in enumerate(wells):
        w = "{}{}".format(w[0], w[1])
        row = int(py.floor(idx / 12))
        col = (idx % 12)
        for c in data:
            shapes = ['--', '*--', '.--']
            for cidx, rep in enumerate(data[c][w]):
                clr = colors[(cidx % 3)]
                shp = shapes[(cidx % 3)]
                curve = data[c][w][rep].rawcurve
                axarr[row, col].plot(time, curve, '{}{}'.format(clr, shp), linewidth=1.5)
                # Only plot tick marks on first column and last row
                if col != 0:
                    axarr[row, col].set_yticks([])
                if row != 7:
                    axarr[row, col].set_xticks([])
                axarr[row, col].tick_params(axis='both', left='on', right='off', bottom='on', top='off')
                axarr[row, col].set_ylim((0, 1.2))
                axarr[row, col].set_title(w)
    f.set_size_inches(35, 18)
    plt.savefig('{}/growthcurves.png'.format(outDir), dpi=100, bbox_inches='tight')




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
parser.add_argument('-p', '--noplate', action='store_true',
                    help='Input wells are not based on a plate')
parser.add_argument('-m', '--heatmap', action='store_true',
                    help='Create growth level heatmap')

args = parser.parse_args()
inputFile = args.infile
outSuffix = args.outsuffix if args.outsuffix else 'out'
outDir = args.outdir
filterFlag = args.filter
newGrowthFlag = args.newgrowth
verbose = args.verbose
debugOut = args.debug
noPlate = args.noplate
hmap = args.heatmap

###############################################################################
# Data Processing
###############################################################################

# Parse data file
printStatus('Parsing input file...')
pmData = PMData.PMData(inputFile, noPlate)
printStatus('Parsing complete.')
if verbose:
    if noPlate:
        printStatus('Plate option not given.')
        printStatus('Found {} samples.'.format(pmData.numClones))
    else:
        printStatus('Found {} samples and {} growth conditions.'.format(
            pmData.numClones, pmData.numConditions))
if debugOut:
    # Print out number of replicates for each clone
    for c, reps in pmData.replicates.items():
        printStatus('DEBUG: {} has {} replicates.'.format(c, len(reps)))


# Perform filter
if filterFlag:
    printStatus('Performing filtering...')
    # Iterate through clones
    for c, repDict in pmData.dataHash.items():
        # Iterate through replicates
        for rep, wellDict in repDict.items():
            # Iterate through wells
            for w, odDict in wellDict.items():
                # Perform filter check
                curveFilter(c, rep, w, odDict['od'], pmData)

    printStatus('Filtering complete.')
    if verbose:
        printStatus('Filtered {} growth curves.'.format(pmData.numFiltered))

elif verbose:
    printStatus('Filtering option not given -- no filtering performed.')


# Create growth curves and logistic models
printStatus('Processing growth curves and creating logistic models...')
finalDataReps = {}
finalDataMean = {}
# Iterate through clones
for c in pmData.clones:
    finalDataReps[c] = {}
    finalDataMean[c] = {}

    # Iterate through media sources
    for w in pmData.wells:
        if not noPlate:
            (ms, gc) = pmData.wells[w]

        finalDataReps[c][w] = {}
        finalDataMean[c][w] = {}

        # Iterate through replicates and determine logistic parameters for each
        tempRepData = py.array([], ndmin=2)
        first = True
        for rep in pmData.replicates[c]:
            if debugOut and noPlate:
                printStatus('DEBUG: Processing {}\t{}\t{}.'.format(c, rep, w))
            elif debugOut:
                printStatus('DEBUG: Processing {}\t{}\t{}\t{}\t{}.'.format(c, rep, ms, gc, w))

            # Create GrowthCurve object for sample
            currCurve = pmData.getODCurve(c, w, rep)
            gCurve = GrowthCurve.GrowthCurve(currCurve, pmData.time)

            # Save sample GrowthCurve object for records
            finalDataReps[c][w][rep] = gCurve

            # Add to temporary multi-dim array for calculating mean
            if first:
                tempRepData = py.array([gCurve.y0, gCurve.asymptote, gCurve.maxGrowthRate,
                                        gCurve.lag, gCurve.growthLevel], ndmin=2)
                first = False
            else:
                add = [gCurve.y0, gCurve.asymptote, gCurve.maxGrowthRate, gCurve.lag, gCurve.growthLevel]
                tempRepData = py.concatenate((tempRepData, [add]))

            if debugOut:
                msg = 'a={}, mgr={}, lag={}'.format(gCurve.asymptote,
                                                    gCurve.maxGrowthRate,
                                                    gCurve.lag)
                printStatus('DEBUG: parameters for {} {} {}: {}'.format(c, rep, w, msg))
        # Create mean logistic curve from replicates' parameters
        meanParams = py.mean(tempRepData, axis=0)
        finalDataMean[c][w]['curve'] = GrowthCurve.logistic(pmData.time, *meanParams[0:4])
        finalDataMean[c][w]['params'] = meanParams
printStatus('Processing complete.')


###############################################################################
# Output Files
###############################################################################

printStatus('Printing output files...')
# Print out filtered data if set
if filterFlag:
    printFiltered(pmData)

# Print out plate info accordingly
if noPlate:
    plateInfo = 'well'
else:
    plateInfo = 'mainsource\tsubstrate\twell'

# logistic_params_sample file: logistic curve parameters for each sample
fhLPSample = open('{}/logistic_params_sample_{}.txt'.format(outDir, outSuffix), 'w')
fhLPSample.write('sample\treplicate\t{}\ty0\tlag\t'.format(plateInfo))
fhLPSample.write('maximumgrowthrate\tasymptote\tgrowthlevel\tsse\n')

# logistic_curves_sample file: logistic curves for each sample
fhLCSample = open('{}/logistic_curves_sample_{}.txt'.format(outDir, outSuffix), 'w')
fhLCSample.write('sample\treplicate\t{}\t'.format(plateInfo))
fhLCSample.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCSample.write('\n')

# logistic_params_mean file. logistic curve parameters (mean)
fhLPMean = open('{}/logistic_params_mean_{}.txt'.format(outDir, outSuffix), 'w')
fhLPMean.write('sample\t{}\ty0\tlag\t'.format(plateInfo))
fhLPMean.write('maximumgrowthrate\tasymptote\tgrowthlevel\n')


# logistic_curves_mean file. logistic curve (mean)
fhLCMean = open('{}/logistic_curves_mean_{}.txt'.format(outDir, outSuffix), 'w')
fhLCMean.write('sample\t{}\t'.format(plateInfo))
fhLCMean.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCMean.write('\n')

# Sort well numbers for print out
if noPlate:
    ws = [(x[0], int(x[1:])) for x in pmData.wells]
else:
    ws = [(x[0], int(x[1:])) for x in pmData.wells.keys()]
sortW = sorted(ws, key=operator.itemgetter(0, 1))
# Iterate through clones
for c, wellDict in finalDataReps.items():
    # Iterate through wells
    for w in sortW:
        w = "{}{}".format(w[0], w[1])
        if noPlate:
            pInfo = w
        else:
            (ms, gc) = pmData.wells[w]
            pInfo = '{}\t{}\t{}'.format(ms, gc, w)

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
        fhLPMean.write('\t'.join(['{:.3f}'.format(x)
                                  for x in (y0, lag, mgr, asymptote, gLevel)]))
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
            fhLPSample.write('\t'.join(['{:.3f}'.format(x)
                                    for x in (y0, lag, mgr, asymptote,
                                              gLevel, sse)]))
            fhLPSample.write('\n')

            # Print logistic curve
            fhLCSample.write('\t'.join(['{:.3f}'.format(x)
                                        for x in curve.dataLogistic]))
            fhLCSample.write('\n')

fhLPSample.close()
fhLCSample.close()
fhLPMean.close()
fhLCMean.close()
printHeatMap(finalDataMean, pmData.clones, sortW, outDir)
curvePlot(finalDataReps, sortW, pmData.time)

printStatus('Printing complete.')
printStatus('Analysis complete.')
