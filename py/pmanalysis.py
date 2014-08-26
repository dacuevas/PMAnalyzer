#!/usr/bin/python
# pmanalysis.py
# Main driver for the phenotype microarray analysis
#
# Author: Daniel A Cuevas
# Created on 22 Nov. 2013
# Updated on 26 Aug. 2014

import argparse
import sys
import time
import datetime
import PMData
import GrowthCurve
import operator
import pylab as py


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

args = parser.parse_args()
inputFile = args.infile
outSuffix = args.outsuffix if args.outsuffix else 'out'
outDir = args.outdir
filterFlag = args.filter
newGrowthFlag = args.newgrowth
verbose = args.verbose
debugOut = args.debug

###############################################################################
# Data Processing
###############################################################################

# Parse data file
printStatus('Parsing input file...')
pmData = PMData.PMData(inputFile)
printStatus('Parsing complete.')
if verbose:
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
    for w, (ms, gc) in pmData.wells.items():
        finalDataReps[c][w] = {}
        finalDataMean[c][w] = {}

        # Iterate through replicates and determine logistic parameters for each
        tempRepData = py.array([], ndmin=2)
        first = True
        for rep in pmData.replicates[c]:
            if debugOut:
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

# logistic_params_sample file: logistic curve parameters for each sample
fhLPSample = open('{}/logistic_params_sample_{}.txt'.format(outDir, outSuffix), 'w')
fhLPSample.write('sample\treplicate\tmainsource\tsubstrate\twell\ty0\tlag\t')
fhLPSample.write('maximumgrowthrate\tasymptote\tgrowthlevel\tsse\n')

# logistic_curves_sample file: logistic curves for each sample
fhLCSample = open('{}/logistic_curves_sample_{}.txt'.format(outDir, outSuffix), 'w')
fhLCSample.write('sample\treplicate\tmainsource\tsubstrate\twell\t')
fhLCSample.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCSample.write('\n')

# logistic_params_mean file. logistic curve parameters (mean)
fhLPMean = open('{}/logistic_params_mean_{}.txt'.format(outDir, outSuffix), 'w')
fhLPMean.write('sample\tmainsource\tsubstrate\twell\ty0\tlag\t')
fhLPMean.write('maximumgrowthrate\tasymptote\tgrowthlevel\n')


# logistic_curves_mean file. logistic curve (mean)
fhLCMean = open('{}/logistic_curves_mean_{}.txt'.format(outDir, outSuffix), 'w')
fhLCMean.write('sample\tmainsource\tsubstrate\twell\t')
fhLCMean.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLCMean.write('\n')

# Sort well numbers for print out
ws = [(x[0], int(x[1:])) for x in pmData.wells.keys()]
sortW = sorted(ws, key=operator.itemgetter(0, 1))
# Iterate through clones
for c, wellDict in finalDataReps.items():
    # Iterate through wells
    for w in sortW:
        w = "{}{}".format(w[0], w[1])
        (ms, gc) = pmData.wells[w]

        # Process mean information
        try:
            curr = finalDataMean[c][w]
        except KeyError:
            # KeyError occurs when all replicates were filtered out so does not
            # exist in final hash
            continue

        # Print sample information
        fhLPMean.write('{}\t{}\t{}\t{}\t'.format(c, ms, gc, w))
        fhLCMean.write('{}\t{}\t{}\t{}\t'.format(c, ms, gc, w))

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
            fhLPSample.write('{}\t{}\t{}\t{}\t{}\t'.format(c, rep, ms, gc, w))
            fhLCSample.write('{}\t{}\t{}\t{}\t{}\t'.format(c, rep, ms, gc, w))

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

printStatus('Printing complete.')
printStatus('Analysis complete.')
