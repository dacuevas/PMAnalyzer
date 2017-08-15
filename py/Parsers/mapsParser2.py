#!/usr/local/bin/python3
# mapsParser2.py
# Parsing script for multi-plate spectrophotometer in ROW format
#
# Author: Daniel A Cuevas
# Created on 24 Aug 2015
# Updated on 14 Aug 2017
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file
# = Each read is a different plate
# = Reads are written in row format (1 row x 96 columns)
# = Each file is a different time point
# === Although, there is a different time point associated with each read
# = First set of rows are column headers (irrelevant)
# = Next row is well information column headers
# === E.g., <TAB> Temperature <TAB> A1 <TAB> A2 <TAB>...
# = Next row is temperature followed by OD reads
# === E.g., <TAB> 27.0 <TAB> 0.07 <TAB> 0.09 <TAB>...
# = After reads is a blank line and an "~End" indicator
# = Finally, there is a text line with the date+timestamp at the end
# === E.g. ...Date Last Saved: 7/8/2015 5:36:31 PMa
# = All of this repeats throughout the file for every plate
# === End of file is reached after every plate is read
# === Number of reads equates to the number of plates
###############################################################

from __future__ import absolute_import, division, print_function
import argparse
import os
import sys
import re
import datetime as dt
from operator import itemgetter

TIMESTAMP = r"(\d+\/\d+\/\d+\s\d+:\d+:\d+\s[AP]M)"


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
    return sortFiles(mypath, files)


def sortFiles(d, files):
    """
    Sort the list of files.
    E.g.: 'Data 1', 'Data 2',..., 'Data 10',...
    """
    # Separate extension from names
    fbases = [os.path.basename(f) for f in files]
    fileNames = [os.path.splitext(f) for f in fbases]
    fileExts = [f[1] for f in fileNames]
    files = [f[0] for f in fileNames]

    # Create a list of tuples of file name and number
    # [("Data", 1), ("Data", 2),..., ("Data", 10),...]
    files2 = [(b[0], int(b[1])) for b in [a.split() for a in files]]
    # Sort by file name then by number
    files2_sorted = sorted(files2, key=itemgetter(0, 1))
    # Recombine and return
    retFiles = []
    for i, fname in enumerate(files2_sorted):
        fn = "{} {}".format(*fname)
        ext = fileExts[i]
        retFiles.append("{}{}{}".format(d, fn, ext))
    return retFiles


def readSampleFile(filepath):
    """
    Parse sample file into dictionary.
    Return sample dictionary and sample order list.
    """
    samples = {}
    sOrder = []
    with open(filepath, "r") as fh:
        for l in fh:
            l = l.rstrip("\n")
            s, r = l.split("\t")
            try:
                samples[s]
            except KeyError:
                samples[s] = []
            samples[s].append(r)
            sOrder.append((s, r))
    return samples, sOrder


def readPlate(mypath):
    """Parse plate file into dictionary"""
    plate = {}
    with open(mypath, "r") as fh:
        for l in fh:
            l = l.rstrip("\n")
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


def readData(data, f, t0, samples, sOrder):
    dataNext = False  # Flag to read in data in next line
    wells = []
    sampleIdx = 0
    with open(f, "r") as fh:
        for ln, l in enumerate(fh):
            # Strip off any blank or newlines on both sides
            l = l.strip()
            if re.match(r"Temperature", l):
                ll = l.split("\t")
                # Check that there are 97 items (1 Temp + 96 wells)
                if len(ll) != 97:
                    errOut("There are not 97 items in header line " + str(ln)
                           + " for file " + f)
                wells = ll[1:]
                dataNext = True

            elif dataNext:
                ll = l.split("\t")
                # Check that there are 97 items (1 Temp + 96 wells)
                if len(ll) != 97:
                    errOut("There are not 97 items in data line " + str(ln))
                odReads = ll[1:]
                dataNext = False

            elif re.search(TIMESTAMP, l):
                # Get sample name and rep
                try:
                    name, rep = sOrder[sampleIdx]
                except ValueError:
                    errOut("Could not extract name and rep in line " + str(ln))
                # Timestamp signifies end of read
                timeObj, t0 = parseTime(l, t0, name, rep)
                data = addOD(data, wells, odReads, timeObj, t0, name, rep)
                sampleIdx += 1


def parseTime(line, t0, name, rep):
    """Parse time from a string and return the datetime object"""
    m = re.search(TIMESTAMP, line)
    timeRaw = m.group(1)
    timeObj = dt.datetime.strptime(timeRaw, "%m/%d/%Y %I:%M:%S %p")

    # Check if this is the first time for the sample or replicate
    if name not in t0:
        t0[name] = {rep: timeObj}
    elif rep not in t0[name]:
        t0[name][rep] = timeObj
    return timeObj, t0


def addOD(data, wells, odReads, timeObj, t0, name, rep):
    """Store data in Data Frame"""
    # Calculate time difference and convert to hours
    tDelta = timeObj - t0[name][rep]
    tDelta = tDelta.days * 24 + tDelta.seconds / 3600
    if tDelta < 0:
        errOut("ERROR: tDelta is negative for " + name + " " + rep)

    for idx, odRead in enumerate(odReads):
        well = parseWell(wells[idx])
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
    # Temporary fix for #ERR in OD
    lastOD = {s: {} for s in data.keys()}

    for s in data:
        # s = sample name
        for r in data[s]:
            if r not in lastOD[s]:
                lastOD[s][r] = {}
            # r = replicate
            for w in sorted(data[s][r], key=itemgetter(0, 1)):
                if w not in lastOD[s][r]:
                    lastOD[s][r][w] = "0.0"
                # w = well tuple: sorted on row then column
                for t, od in sorted(data[s][r][w].items()):
                    if od == "#ERR":
                        od = lastOD[s][r][w]
                    lastOD[s][r][w] = od
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
parser.add_argument("sample", help="Sample and replicate file listing")
parser.add_argument("indir", help="Directory containing data files")
parser.add_argument("-p", "--plate", help="Plate file for wells")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Increase output for status messages")

args = parser.parse_args()
inDir = args.indir
sampleFile = args.sample
plateFile = args.plate
verbose = args.verbose


###############################################################################
# Preprocessing files
###############################################################################

# Get data files
filenames = getDataFiles(inDir)

# Get sample and replicate list
samples, sOrder = readSampleFile(sampleFile)

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
    readData(data, f, t0, samples, sOrder)

# Print out header line
if plateFile:
    print("sample\trep\tmainsource\tcompound\twell\tplate\ttime\tod")
else:
    print("sample\trep\twell\ttime\tod")
printData(data, plate, plateName)
