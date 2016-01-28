#!/usr/loca/bin/python3
# info_vs_time.py
# Parsing script for files that include all information with time points
# as column headers.
#
# Author: Daniel A Cuevas
# Created on 23 Jul 2015
# Updated on 20 Aug 2015
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file, spreadsheet format
# = First row are column headers
# = First set of columns are information
# === E.g., sample   mainsource   substrate   plate   well
# === Order will not matter using Pandas library
#
# = After information columns come time points
# === E.g., plate   well   0   0.5   1.0   1.5
#
# = Proceeding rows are data for each column
###############################################################

from __future__ import absolute_import, division, print_function
import argparse
import os
import sys
import pandas as pd


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


def readData(data, f):
    """Use Pandas library to melt data file into appropriate format"""
    d = pd.read_csv(f, delimiter="\t")
    id_vars = ["sample", "rep", "mainsource", "compound", "well", "plate"]
    var_name = "time"
    value_name = "od"
    dMelt = pd.melt(d, id_vars=id_vars, var_name=var_name,
                    value_name=value_name)
    # Drop rows with NaN
    # This happens when data is not recorded for a sample at a timepoint
    # when another sample has a data
    dMelt.dropna(inplace=True)

    dMelt.set_index(id_vars, inplace=True)
    if data is None:
        return dMelt
    else:
        return pd.concat([data, dMelt])


###############################################################################
# Argument parsing
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument("indir", help="Directory containing data files")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="Increase output for status messages")

args = parser.parse_args()
inDir = args.indir
verbose = args.verbose

# Get data files
filenames = getDataFiles(inDir)

###############################################################################
# Data extraction files
###############################################################################

# All data will be contained in a dictionary
# Keys will be ordered as follows:
#     1) Sample name - str
#     2) Replicate - str
#     3) Well - 2-tuple (str, int)
#     4) Plate - str
#     5) Time - float
#     6) OD value - str
data = None

# Read in data for each file
for f in filenames:
    data = readData(data, f)

# Print out results
print(data.to_csv(None, sep="\t", header=True, index=True), end="")
