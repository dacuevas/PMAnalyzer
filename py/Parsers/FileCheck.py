#!/usr/local/bin/python3
# FileCheck.py
# Perform a file format check for the PMAnalyzer pipeline
#
# Author: Daniel A Cuevas
# Created on 27 Dec 2016
# Updated on 09 Nov 2017

from __future__ import absolute_import, print_function
import argparse
import os
import sys
import re


###############################################################################
# Utility methods
###############################################################################
def errOut(msg):
    """Send a system exit exception with the given message"""
    script = os.path.basename(__file__)  # No longer using script name
    sys.exit("ERROR during file format check.\n    {}".format(msg))


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


def getSampleAndRep(file):
    """Extract sample name and replicate ID from file names"""
    m = re.match(r"^([A-Za-z0-9-.]+)_([A-Za-z0-9]+)", file)
    if m is None:
        errOut("Improper filename for parser. Could not extract name and replicate from filename: "
               "{}".format(os.path.basename(file)))
    name, rep = m.group(1, 2)
    return name, rep


def checkFormat1(file):
    """Check formatting of text in file format 1: mapsParser1.py"""
    tsre = r"(\d+\/\d+\/\d+\s\d+:\d+:\d+\s[AP]M)\s+"
    wellre = r"^\w\d+\s+[0-9.]+$"
    with open(file) as f:
        foundTimestamp = False  # Flag when timestamp found
        foundData = False  # Flag when OD reads found
        numData = 0
        for lnum, l in enumerate(f, start=1):
            l = l.rstrip("\n")
            # Start of file should be a timestamp in the format of:
            # MM/DD/YYY HH:MM:SS [AM|PM]
            if lnum == 1:
                m = re.match(tsre, l)
                if m == None:
                    errOut("Start of file not a timestamp for filename: "
                           "{}".format(os.path.basename(file)))
                foundTimestamp = True

            # Looking for next timestamp
            elif not foundTimestamp and not foundData and re.match(tsre, l):
                foundTimestamp = True
                numData += 1

            # Looking if next timestamp found but data was never found
            elif foundTimestamp and not foundData and re.match(tsre, l):
                msg = "File {} had incorrect format.".format(file)
                msg += " Error encountered at line {}.".format(lnum)
                msg += " No data found between timestamps."
                msg += " Data may be in incorrect format."
                errOut(msg)

            # Looking for OD data after timestamp
            elif foundTimestamp and re.match(wellre, l):
                foundData = True

            # Looking for end of OD data
            elif foundTimestamp and foundData and not re.match(wellre, l):
                foundTimestamp = False
                foundData = False

        numData += 1

        # Should be more than 2 data points
        if numData < 3 and lnum > 600:
            errOut("Less than 3 data points found but file is longer than 600 "
                   "lines long. File may be incorrectly formatted.")


def checkFormat2(file):
    """Check formatting of text in file format 2: info_vs_time.py"""
    with open(file) as f:
        numColumns = 0
        for lnum, l in enumerate(f, start=1):
            l = l.rstrip("\n")
            ll = l.split("\t")
            if lnum == 1:
                numColumns = len(ll)
                # Check number of columns, must be at least 7
                if numColumns < 7:
                    errOut("Missing columns in file. At least 7 are required. "
                           "Please view help section for correct info.")

                # Check column info headers
                if (ll[0] != "sample" or ll[1] != "rep" or
                    ll[2] != "mainsource" or ll[3] != "compound" or
                    ll[4] != "plate" or ll[5] != "well"):
                    errOut("Column names are incorrect. "
                           "Please view help section for correct info.")

                # Check time, must be only numbers and decimal points
                for t in ll[6:]:
                    if not re.match(r"^[0-9.]+$", t):
                        errOut("Only numbers allowed in time headers. "
                               "See time header '{}'.".format(t))
                continue

            # All other data points
            # Check number of columns
            if len(ll) != numColumns:
                errOut("Inconsistent number of columns in data file "
                       "line {}.".format(lnum))

            # Check data points are not empty
            for data in ll:
                if data == "":
                    errOut("Missing data field found in data file "
                           "line {}".format(lnum))

            # Check OD data points are only floating point numbers
            for data in ll[6:]:
                if not re.match(r"^[0-9.]+$", data):
                    errOut("OD data must be floating point numbers only. See "
                           "line {}".format(lnum))


def checkFormat3(file):
    """Check formatting of text in file format 3: mapsParser2.py"""
    # Implementation coming soon
    pass


def checkFormat4(file):
    """Check formatting of text in file format 4: mapsParser3.py"""
    # Implementation coming soon
    pass


def checkFormat5(file):
    """Check formatting of text in file format 5: mapsParser4.py"""
    # Implementation coming soon
    pass


def checkFormat6(file):
    """Check formatting of text in file format 6: mapsParser5.py"""
    # Implementation coming soon
    pass


###############################################################################
# Argument parsing
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument("indir", help="Directory containing data files")
parser.add_argument("file_type", help="Number for the file parser to use",
                    type=int, choices=range(1, 7))
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Increase output for status messages")

args = parser.parse_args()
inDir = args.indir
fileType = args.file_type
verbose = args.verbose

###############################################################################
# Processing files
###############################################################################
# Get data files
filenames = getDataFiles(inDir)

for f in filenames:
    # Attempt to extract sample names and replicate IDs
    if fileType == 1:
        name, rep = getSampleAndRep(os.path.basename(f))

    # Check format of file
    if fileType == 1:
        checkFormat1(f)
    elif fileType == 2:
        checkFormat2(f)
    elif fileType == 3:
        checkFormat3(f)
    elif fileType == 4:
        checkFormat4(f)
    elif fileType == 5:
        checkFormat5(f)
    elif fileType == 6:
        checkFormat6(f)
