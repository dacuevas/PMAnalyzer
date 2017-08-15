#!/usr/local/bin/python3
# mapsParser1.py
# Parsing script for GT Analyst multi-plate spectrophotometer
#
# Author: Daniel A Cuevas
# Created on 13 Apr 2015
# Updated on 14 Aug 2017

from __future__ import absolute_import, division, print_function
import argparse
import os
import sys
import re
import datetime as dt
from operator import itemgetter


###############################################################################
# Utility methods
###############################################################################
def errOut(msg, parser=None):
    """Send a system exit exception with the given message"""
    if parser is not None:
        parser.print_help()
    script = os.path.basename(__file__)
    sys.exit("ERROR in {}\n    {}".format(script, msg))


def getDataFiles(mypath):
    """Check directory exists and return list of non-hidden files"""
    # Add forward slash if directory does not have one
    if not mypath.endswith("/"):
        mypath = mypath + "/"
    # Check if directory exists
    if not os.path.isdir(mypath):
        errOut("Data directory not found: {}".format(mypath))
    files = ["{}{}".format(mypath, f) for f in os.listdir(mypath)
             if os.path.isfile(os.path.join(mypath, f)) and
             not f.startswith(".")]
    return files


def readPlate(mypath):
    """Parse plate file into dictionary"""
    plate = {}
    with open(mypath, "r") as fh:
        for l in fh:
            l = l.rstrip()
            # Order of elements are:
            #   1) well
            #   2) mainsource
            #   3) compound
            #   4) concentration
            well, ms, cpd, conc = l.split("\t")
            well = parseWell(well)
            plate[well] = (ms, cpd)
    return plate


def parseWell(well):
    """Convert well string to 2-tuple"""
    well = re.match(r"(\w)(\d+)", well).group(1, 2)
    well = (well[0], int(well[1]))
    return well


def readData(data, f, t0):
    # Get only filename
    fname = os.path.basename(f)
    # Remove file extension
    fname = os.path.splitext(fname)[0]

    # Extract sample name and replicate from file
    m = re.match(r"^([A-Za-z0-9-.]+)_([A-Za-z0-9]+)", fname)
    if m is None:
        errOut("Could not extract name and replicate from filename: "
               "{}".format(f))
    name, rep = m.group(1, 2)

    # Be sure to ignore the BACKGROUND data points in the file
    inBackground = False
    with open(f, "r") as fh:
        for l in fh:
            l = l.rstrip()
            if re.match(r"BACKGROUND", l):
                # Check for background data - this always comes after raw data
                inBackground = True

            elif re.match(r"\d.*?\s.{11}\s+", l):
                # No longer in background data
                inBackground = False
                timeObj, t0 = parseTime(l, t0, name, rep)

            elif not inBackground and re.match(r"\w+\s+[0-9.]+", l):
                # Extract well and data
                data = parseOD(l, data, t0, timeObj, name, rep)


def parseTime(line, t0, name, rep):
    """Parse time from a string and return the datetime object"""
    m = re.match(r"(\d+\/\d+\/\d+\s\d+:\d+:\d+\s[AP]M)\s+", line)
    timeRaw = m.group(1)
    timeObj = dt.datetime.strptime(timeRaw, "%m/%d/%Y %I:%M:%S %p")

    # Check if this is the first time for the sample or replicate
    if name not in t0 or rep not in t0[name]:
        # Store starting time
        t0[name] = {rep: timeObj}
    return timeObj, t0


def parseOD(line, data, t0, timeObj, name, rep):
    """Parse data from a string and store in Data Frame"""
    # Extract well and data
    m = re.match(r"(\w+)\s+([0-9.]+)", line)
    well, odRead = m.group(1, 2)
    well = parseWell(well)

    # Calculate time difference and convert to hours
    tDelta = timeObj - t0[name][rep]
    tDelta = tDelta.days * 24 + tDelta.seconds / 3600

    # Check if data keys are initialized
    if name not in data:
        data[name] = {rep: {well: {tDelta: odRead}}}
    elif rep not in data[name]:
        data[name][rep] = {well: {tDelta: odRead}}
    elif well not in data[name][rep]:
        data[name][rep][well] = {tDelta: odRead}
    else:
        data[name][rep][well][tDelta] = odRead
    return data


def printData(data, plate=None, pn=None):
    """Print parsed data from dictionary. Output is tab-delimited.
    Columns:
        1) sample
        2) replicate
        3) main source*
        4) compound*
        5) well#
        6) plate name*
        7) time
        8) OD reading
    *Only printed when given a plate file
    """
    for s in data:
        # s = sample name
        for r in data[s]:
            # r = replicate
            for w in sorted(data[s][r], key=itemgetter(0, 1)):
                # w = well tuple: sorted on row then column
                for t, od in sorted(data[s][r][w].items()):
                    # t = time
                    wId = "".join((w[0], str(w[1])))
                    if plate:
                        ms, cpd = plate[w]
                        print("\t".join((s, r, ms, cpd, wId, pn,
                                         "{:.1f}".format(t), str(od))))
                    else:
                        print("\t".join((s, r, wId,
                                         "{:.1f}".format(t), str(od))))


###############################################################################
# Argument parsing
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument("indir", help="Directory containing data files")
parser.add_argument("-o", "--outsuffix",
                    help="Suffix appended to output files. Default is 'out'")
parser.add_argument("-p", "--plate", help="Plate file for wells")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Increase output for status messages")

args = parser.parse_args()
inDir = args.indir
outSuffix = args.outsuffix if args.outsuffix else "out"
verbose = args.verbose
plateFile = args.plate


###############################################################################
# Preprocessing files
###############################################################################

# Get data files
filenames = getDataFiles(inDir)

# Check that the plate file exists
# If not, check if it exists in PMAnalyzer premade plates
if plateFile and not os.path.isfile(plateFile):
    d = os.path.expanduser("~") + "/Projects/PMAnalyzer/plates/"
    plateFile = d + plateFile
    if not os.path.isfile(plateFile):
        errOut("Plate file not found: {}".format(plateFile))

if plateFile:
    plate = readPlate(plateFile)
    # Get only filename
    plateName = os.path.basename(plateFile)
    # Remove file extension
    plateName = os.path.splitext(plateName)[0]
else:
    plate = None
    plateName = None

###############################################################################
# Data extraction files
###############################################################################

# All data will be contained in a dictionary
# Keys will be ordered as follows:
#     1) Sample name - str
#     2) Replicate - str
#     3) Well - 2-tuple (str, int)
#     4) Time - float
#     5) OD value - str
data = {}
t0 = {}  # Dictionary containing starting time points for each replicate

# Read in data for each file
for f in filenames:
    readData(data, f, t0)

# Print out header line
if plateFile:
    print("sample\trep\tmainsource\tcompound\twell\tplate\ttime\tod")
else:
    print("sample\trep\twell\ttime\tod")
printData(data, plate, plateName)
