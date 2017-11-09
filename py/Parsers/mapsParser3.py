#!/usr/local/bin/python3
# mapsParser3.py
# Parsing script for multi-plate spectrophotometer in PLATE format
#
# Author: Daniel A Cuevas
# Created on 08 Mar 2016
# Updated on 09 Mar 2016
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file
# = Each read is a different plate
# = Reads are written in plate format (8 rows x 12 columns)
# = Each file is a different plate cycle
# === Although, there is a different time point associated with each read
# = First set of rows are column headers (irrelevant)
# = Next row is the 12 well information column headers
# === E.g., <TAB> Temperature <TAB> 1 <TAB> 2 <TAB>...
# = Next 8 rows is temperature followed by OD reads for row A through H
# === E.g., <TAB> 27.0 <TAB> 0.07 <TAB> 0.09 <TAB>...
# ===       <TAB> Blank <TAB> 0.09 <TAB> 0.08 <TAB>...
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
from itertools import product

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
    data_row = 9  # Will only read data for data row 1 through 8
    sampleIdx = 0
    od_reads = []
    with open(f, "r") as fh:
        for ln, l in enumerate(fh, start=1):
            # Strip off any blank or newlines on both sides
            l = l.strip()
            if re.match(r"Temperature", l):
                ll = l.split("\t")
                # Check that there are 13 items (1 Temp + 12 columns)
                if len(ll) != 13:
                    errOut("There are not 13 items in header line " + str(ln)
                           + " for file " + f)
                data_row = 1

            elif data_row <= 8:
                ll = l.split("\t")
                # Check that there are 12 items (12 columns)
                # First data row has 13 items -- first item is the temperature
                if data_row == 1:
                    if len(ll) != 13:
                        errOut("There are not 13 items in the first data line "
                               + str(ln) + " for file " + f)
                    od_reads.append(ll[1:])
                elif data_row > 1:
                    if len(ll) != 12:
                        errOut("There are not 12 items in data line " + str(ln)
                               + " for file " + f)
                    od_reads.append(ll)
                data_row += 1

            elif re.search(TIMESTAMP, l):
                # Get sample name and rep
                try:
                    name, rep = sOrder[sampleIdx]
                except ValueError:
                    errOut("Could not extract name and rep in line " + str(ln))
                # Timestamp signifies end of read
                timeObj, t0 = parseTime(l, t0, name, rep)
                data = addOD(data, od_reads, timeObj, t0, name, rep)
                sampleIdx += 1
                od_reads = []


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


def addOD(data, od_reads, timeObj, t0, name, rep):
    """Store data in Data Frame"""
    # Calculate time difference and convert to hours
    tDelta = timeObj - t0[name][rep]
    tDelta = tDelta.days * 24 + tDelta.seconds / 3600
    if tDelta < 0:
        errOut("ERROR: tDelta is negative for " + name + " " + rep)

    wells = ["{}{}".format(*w) for w in product("ABCDEFGH", range(1, 13))]

    # Iterate through reads to add
    # First check if od_reads is the correct number of rows
    if len(od_reads) != 8:
        errOut("ERROR: od_reads does not contain 8 rows of data. There are"
               " only " + str(len(od_reads)) + " rows. Name: " + name
               + " Rep: " + rep)
    for row_idx, od_row in enumerate(od_reads):
        # Check that the row contains 12 items of data
        if len(od_row) != 12:
            errOut("ERROR: od_row " + str(row_idx) + " does not contain "
                   "12 columns of data. Only " + str(len(od_reads)) + " cols."
                   " Name: " + name + " Rep: " + rep)
        for col_idx, od_read in enumerate(od_row):
            idx = row_idx * 12 + col_idx
            well = parseWell(wells[idx])
            # Check if data keys are initialized
            if name not in data:
                data[name] = {rep: {well: {tDelta: od_read}}}
            elif rep not in data[name]:
                data[name][rep] = {well: {tDelta: od_read}}
            elif well not in data[name][rep]:
                data[name][rep][well] = {tDelta: od_read}
            else:
                data[name][rep][well][tDelta] = od_read
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