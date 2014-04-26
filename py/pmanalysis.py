#!/usr/bin/python
# pmanalysis.py
# Main driver for the phenotype microarray analysis
#
# Author: Daniel A Cuevas
# Created on 22 Nov. 2013
# Updated on 10 Apr. 2014

import argparse
import sys
import time
import datetime
import PMData
import GrowthCurve
import operator


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
logData = {}
# Iterate through clones
for c in pmData.clones:
    logData[c] = {}

    # Iterate through media sources
    for w, (ms, gc) in pmData.wells.items():
        if debugOut:
            printStatus('DEBUG: Processing {}\t{}\t{}\t{}.'.format(c, ms, gc, w))
        logData[c][w] = {}
        curves = pmData.getCloneReplicates(c, w, filterFlag)

        # Add curve to logData hash
        # Will not add if:
        # 1. Filtering is on
        # 2. All replicates were filtered out
        if len(curves) > 0:
            gc = GrowthCurve.GrowthCurve(curves, pmData.time)
            logData[c][w] = gc
            if debugOut:
                msg = 'a={}, mgr={}, lag={}'.format(gc.asymptote,
                                                    gc.maxGrowthRate,
                                                    gc.lag)
                printStatus('DEBUG: parameters for {} {}: {}'.format(c, w, msg))
printStatus('Processing complete.')


###############################################################################
# Output Files
###############################################################################

printStatus('Printing output files...')
# Print out filtered data if set
if filterFlag:
    printFiltered(pmData)

# curveinfo file: curve parameters for each sample
fhInfo = open('{}/curveinfo_{}.txt'.format(outDir, outSuffix), 'w')
fhInfo.write('sample\tmainsource\tsubstrate\twell\tlag\t')
fhInfo.write('maximumgrowthrate\tasymptote\tgrowthlevel\n')

# logistic_curve file: logistic curves
fhLogCurve = open('{}/logistic_curves_{}.txt'.format(outDir, outSuffix), 'w')
fhLogCurve.write('sample\tmainsource\tsubstrate\twell\t')
fhLogCurve.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhLogCurve.write('\n')

# median file: median curves
fhMedCurve = open('{}/median_curves_{}.txt'.format(outDir, outSuffix), 'w')
fhMedCurve.write('sample\tmainsource\tsubstrate\twell\t')
fhMedCurve.write('\t'.join(['{:.1f}'.format(x) for x in pmData.time]))
fhMedCurve.write('\n')

# Sort well numbers for print out
ws = [(x[0], int(x[1:])) for x in pmData.wells.keys()]
sortW = sorted(ws, key=operator.itemgetter(0, 1))
# Iterate through clones
for c, wellDict in logData.items():
    # Iterate through wells
    for w in sortW:
        w = "{}{}".format(w[0], w[1])
        try:
            curve = wellDict[w]
        except KeyError:
            # KeyError occurs when all replicates were filtered out so does not
            # exist in final hash
            continue
        (ms, gc) = pmData.wells[w]

        # Print sample information
        fhInfo.write('{}\t{}\t{}\t{}\t'.format(c, ms, gc, w))
        fhLogCurve.write('{}\t{}\t{}\t{}\t'.format(c, ms, gc, w))
        fhMedCurve.write('{}\t{}\t{}\t{}\t'.format(c, ms, gc, w))

        # Print OD readings
        lag = curve.lag
        mgr = curve.maxGrowthRate
        asymptote = curve.asymptote
        gLevel = curve.growthLevel
        fhInfo.write('\t'.join(['{:.3f}'.format(x)
                                for x in (lag, mgr, asymptote, gLevel)]))
        fhInfo.write('\n')

        # Print logistic curves
        fhLogCurve.write('\t'.join(['{:.3f}'.format(x)
                                    for x in curve.dataLogistic]))
        fhLogCurve.write('\n')

        # Print out median curve
        fhMedCurve.write('\t'.join(['{:.3f}'.format(x)
                                    for x in curve.dataMed]))
        fhMedCurve.write('\n')

fhInfo.close()
fhLogCurve.close()
fhMedCurve.close()

printStatus('Printing complete.')
printStatus('Analysis complete.')
