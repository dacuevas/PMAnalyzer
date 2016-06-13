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
# First check if known file name is given or use full path
if not os.path.isfile(args.plate):
    print("Checking PMAnalyzer directory for plate...", file=sys.stderr)
    d = os.path.expanduser("~") + "/Projects/PMAnalyzer/plates/"
    args.plate = d + args.plate
with open(args.plate, "r") as f:
    for l in f:
        w, ms, c, conc = l.strip().split("\t")
        plateinfo[w] = (ms, c)

with open(args.input, "r") as f:
    header = next(f)
    header = header.strip()
    ll = header.split("\t")
    wellIdx = ll.index("well")
    if args.nocolumn:
        # Plate information does not exist yet
        # Print out new columns
        header = "\t".join(ll[0:2])
        header += "\tmainsource\tcompound\tplate\t"
        header += "\t".join(ll[2:])
    print(header)
    for l in f:
        l = l.strip()
        ll = l.split("\t")
        w = ll[wellIdx]
        if args.nocolumn:
            # Plate information does not exist yet
            # Add new columns
            p = re.search(r"([A-z0-9-_.]+).txt$", args.plate).group(1)
            line = "\t".join(ll[0:2]) + "\t"
            line += "\t".join([ms, c, p]) + "\t"
            line += "\t".join(ll[2:])
        else:
            # Plate information already exists
            # Replace mainsource and compound
            ms, c = plateinfo[w]
            line = "\t".join(ll[0:2]) + "\t"
            line += "\t".join([ms, c]) + "\t"
            line += "\t".join(ll[4:])
        print(line)

print("Complete!", file=sys.stderr)
