#!/usr/bin/python2.7
from __future__ import print_function
import argparse
import re
import sys
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input PMAnalyzer file")
parser.add_argument("plate", help="Plate file name")
parser.add_argument("-n", "--nocolumn", action="store_true",
                    help="Plate information columns do not exist")

args = parser.parse_args()

plateinfo = {}
plate = re.search(r"([A-z0-9-_.]+).txt$", args.plate).group(1)
# First check if known file name is given or use full path
if not os.path.isfile(args.plate):
    print("Checking PMAnalyzer directory for plate...", file=sys.stderr)
    d = os.path.expanduser("~") + "/Projects/PMAnalyzer/plates/"
    args.plate = d + args.plate

# Read in new plate information
with open(args.plate, "r") as f:
    for l in f:
        w, ms, c, conc = l.strip().split("\t")
        plateinfo[w] = (ms, c)

# Edit raw curves file with new plate information
with open(args.input, "r") as f:
    # Read in header
    header = next(f)
    header = header.strip()
    ll = header.split("\t")
    # Find index to well column
    wellIdx = ll.index("well")
    otherColsIdx = 6  # Starting index to other columns (e.g., time, od)
                      # Will be reset below if there is no plate info
    if args.nocolumn:
        # Plate information does not exist yet
        # Print out new columns
        header = "\t".join(ll[0:2])
        header += "\tmainsource\tcompound\tplate\t"
        header += "\t".join(ll[2:])
        otherColsIdx = 3  # No ms, c, and plate info
    print(header)
    for l in f:
        l = l.strip()
        ll = l.split("\t")

        # Grab well and plate info
        w = ll[wellIdx]
        ms, c = plateinfo[w]

        # Print out sample and replicate
        line = "\t".join(ll[0:2]) + "\t"
        # Print out plate info
        line += "\t".join([ms, c, w, plate]) + "\t"
        # Print out rest of columns
        line += "\t".join(ll[otherColsIdx:])
        print(line)

print("Complete!", file=sys.stderr)
