#!/usr/local/bin/python3
# mapsParser4.py
# Parsing script for multi-plate spectrophotometer in ROW format
#
# Author: Daniel A Cuevas
# Created on 11 Nov 2017
# Updated on 11 Nov 2017
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file, spreadsheet format
# = File is for a single plate
# = First two rows are machine-generated information
# = Third row are column headers, a total of 98 columns with 98 tabs
# === Time <TAB> Temperature <TAB> A1 <TAB> A2 <TAB> ... <H12> <TAB> <NEWLINE>
# = Remaining rows are data values for columns
# = Time format: D.HH:MM:SS
# === Regex for time format is (\d\.)?\d{2}:\d{2}:\d{2}
# === Days and period are optional
# = After reads is a blank line and an "~End" indicator
# = Finally, there is a text line with the date+timestamp at the end (not used)
# === E.g. ...Date Last Saved: 7/8/2015 5:36:31 PM
###############################################################

from __future__ import absolute_import, division, print_function
import argparse
import os
import sys
import re
import datetime as dt
from operator import itemgetter
from itertools import product

TIMESTAMP = re.compile(r"((\d)\.)?(\d{2}:\d{2}:\d{2})")


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

    data_rows = False  # Flag when to read data
    with open(f, "r") as fh:
        for ln, l in enumerate(fh, start=1):
            # Strip off any blank or newlines on both sides
            l = l.strip()

            # Check if we reached the end of the file
            if data_rows and len(l) == 0:
                break

            elif re.match(r"Time", l):
                ll = l.split("\t")
                # Check that there are 98 items (1 Time, 1 Temp, 96 wells)
                if len(ll) != 98:
                    errOut("There are not 98 items in header line " + str(ln)
                           + " for file " + f)
                data_rows = True

            elif data_rows:
                ll = l.split("\t")
                # Check that there are 98 items
                if len(ll) != 98:
                    errOut("There are not 98 items in data line " + str(ln)
                           + " for file " + f)
                timeObj, t0 = parseTime(ll[0], t0, name, rep)
                data = addOD(ll[2:], data, timeObj, t0, name, rep)



def parseTime(time_str, t0, name, rep):
    """Parse time from a string and return the datetime object"""
    m = TIMESTAMP.match(time_str)
    #timeRaw = ":".join(list(m.group(2, 3, 4)))
    timeRaw = m.group(3)
    timeObj = dt.datetime.strptime(timeRaw, "%H:%M:%S")
    # Check if there was a day present
    if m.group(2):
        timeObj += dt.timedelta(days=int(m.group(2)))

    # Check if this is the first time for the sample or replicate
    if name not in t0:
        t0[name] = {rep: timeObj}
    elif rep not in t0[name]:
        t0[name][rep] = timeObj
    return timeObj, t0


def addOD(od_reads, data, timeObj, t0, name, rep):
    """Store data in Data Frame"""
    # Calculate time difference and convert to hours
    tDelta = timeObj - t0[name][rep]
    tDelta = tDelta.days * 24 + tDelta.seconds / 3600
    if tDelta < 0:
        errOut("ERROR: tDelta is negative for " + name + " " + rep)

    wells = ["{}{}".format(*w) for w in product("ABCDEFGH", range(1, 13))]

    # Iterate through reads to add
    for idx, od_read in enumerate(od_reads):
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
parser.add_argument("-p", "--plate", help="Plate file for wells")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Increase output for status messages")

args = parser.parse_args()
inDir = args.indir
plateFile = args.plate
verbose = args.verbose


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