#!/usr/bin/python
# pmanalyzer_statistics.py
# Produce statistics from the PMAnalyzer results (logistic_params_sample* file)
#
# Author: Daniel A Cuevas
# Created on 23 Sep 2016
# Updated on 23 Sep 2016

import argparse
import PMUtil as util
import pylab as py
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("results_file",
                    help="logistic_params_sample* results file")
parser.add_argument('outdir',
                    help='Directory to store output files')
parser.add_argument('-o', '--outsuffix',
                    help='Suffix appended to output files. Default is "out"')
parser.add_argument('-p', '--plate', action='store_true',
                    help='Input wells are based on a plate')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Increase output for status messages')
args = parser.parse_args()

fn = args.results_file
outSuffix = args.outsuffix if args.outsuffix else 'out'
outDir = args.outdir
verbose = args.verbose
plateFlag = args.plate
data = pd.read_table(fn, header=0)

if verbose:
    util.printStatus("Reading in {} to output statistics".format(fn))

grp = ["sample", "rep"]
cols = ["y0", "lag", "maxgrowth", "asymptote", "growthlevel", "glscaled",
        "r", "auc_raw", "auc_rshift", "auc_log", "auc_lshift", "err_mse"]
all_average = data.groupby(grp)[cols].mean()
all_median = data.groupby(grp)[cols].median()
all_stdev = data.groupby(grp)[cols].std()

# Rename columns to merge
all_average.columns = [x + "_avg" for x in all_average.columns]
all_median.columns = [x + "_median" for x in all_median.columns]
all_stdev.columns = [x + "_stdev" for x in all_stdev.columns]

# Join DataFrames into one
all_data = all_average.join(all_median.join(all_stdev))

###############################################################################
# Media
###############################################################################
well_grp = ["well"]
if plateFlag:
    well_grp += ["mainsource", "compound"]

well_average = data.groupby(well_grp)[cols].mean()
well_median = data.groupby(well_grp)[cols].median()
well_stdev = data.groupby(well_grp)[cols].std()

# Rename columns to merge
well_average.columns = [x + "_avg" for x in well_average.columns]
well_median.columns = [x + "_median" for x in well_median.columns]
well_stdev.columns = [x + "_stdev" for x in well_stdev.columns]

# Join DataFrames into one
well_data = well_average.join(well_median.join(well_stdev))

###############################################################################
# Output
###############################################################################
col = ["y0_avg", "y0_median", "y0_stdev",
       "lag_avg", "lag_median", "lag_stdev",
       "maxgrowth_avg", "maxgrowth_median", "maxgrowth_stdev",
       "asymptote_avg", "asymptote_median", "asymptote_stdev",
       "growthlevel_avg", "growthlevel_median", "growthlevel_stdev",
       "glscaled_avg", "glscaled_median", "glscaled_stdev",
       "r_avg", "r_median", "r_stdev",
       "auc_raw_avg", "auc_raw_median", "auc_raw_stdev",
       "auc_rshift_avg", "auc_rshift_median", "auc_rshift_stdev",
       "auc_log_avg", "auc_log_median", "auc_log_stdev",
       "auc_lshift_avg", "auc_lshift_median", "auc_lshift_stdev",
       "err_mse_avg", "err_mse_median", "err_mse_stdev"]
if verbose:
    util.printStatus('Saving statistics files.')
all_data.to_csv(
    "{}/sample_statistics_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)
well_data.to_csv(
    "{}/well_statistics_{}.txt".format(outDir, outSuffix),
    sep="\t", header=True, index=True, float_format="%.3f", columns=col)
if verbose:
    util.printStatus('Statistics analysis complete.')
