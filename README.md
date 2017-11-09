## *What is PMAnalyzer?*
PMAnalyzer is a command-line software pipeline that performs bacterial growth curve analysis for **M**ulti-phenotype **A**ssay **P**lates (MAPs).
 MAPs are 96-well micro-titer plates used to measure biomass formation through a spectrophotometer. Through the use of Python's SciPy and NumPy packages,
 PMAnalyzer performs curve-fitting in order to assess lag phase, growth rate, biomass yield, and growth level on each growth curve.

## *Dependencies*
#### Languages
1. Bash
2. R
3. Python 2.7.X+ or 3.X.X

#### R Modules/Packages
1. ggplot2
2. ggthemes
3. reshape2
4. grid
5. plyr

#### Python Modules/Packages
1. SciPy
2. NumPy
3. Matplotlib (for graphs and figures)
4. Pandas

## *Input files*
The parsing scripts from PMAnalyzer are generated from the Molecular Devices Analyst GT multi-plate plate reader. Sample input files can be viewed in the
 `sample/data_csedlakii` directory.

Filenames must be in the format of `[sample]_[replicate]_[other].txt`. For example

| Sample | Replicate | Filename |
|:------:|:---------:|:--------:|
| csedlakii | rep1 | *csedlakii_rep1_2014Nov12.txt* |
| csedlakii | rep2 | *csedlakii_rep2_2014Nov12.txt* |
| csedlakii | rep3 | *csedlakii_rep3_2014Nov12.txt* |
| ecoli | rep1 | *ecoli_rep1_2014Nov15.txt* |
| ecoli | rep2 | *ecoli_rep2_2014Nov15.txt* |
| ecoli | rep3 | *ecoli_rep3_2014Nov15.txt* |

**sample** and **replicate** names can have letters `A-z`, numbers `0-9`, hyphens `-`, and periods `.` (periods for sample only).

## *Output files*
(*All output text files are tab-delimited*)
- Various growth curves text files (raw, mean, median, logistic)
- Growth curve parameter text files (raw, logistic)
- Various growth curve figures
- Growth level heat-map

#### Sample Images
![heatmap](https://github.com/dacuevas/PMAnalyzer/blob/develop/sample/sample_results/growthlevels.png "Growth Level Heatmap")

![growthcurves](https://github.com/dacuevas/PMAnalyzer/blob/develop/sample/sample_results/avg_R.S.3.png "C. sedlakii Growth Curves")

## *Usage*
```
        ___          _               _
       / _ \/\/\    /_\  _ __   __ _| |_   _ _______ _ __
      / /_)/    \  //_\\| '_ \ / _` | | | | |_  / _ \ '__|
     / ___/ /\/\ \/  _  \ | | | (_| | | |_| |/ /  __/ |
     \/   \/    \/\_/ \_/_| |_|\__,_|_|\__, /___\___|_|
                                        |___/

PMAnalyzer version 2.4

usage: runPM -c config_filename [Options]
OR
usage: runPM -i PM_data_directory -d output_directory -o output_name [Options]

Required
   -c [config_filename]    : Configuration file listing all parameters (this
takes precedence)
   OR
   -i [PM_data_directory]  : PM data directory path to files
   -d [output_directory]   : Output directory (will create if non-existent)
   -o [output_name]        : Suffix appended to all output files

Optional
   --clear                 : Clear out results directory of old files
   --debug                 : Print out debug messages
   -g                      : Flag to use experimental growth level calculation. Range: [1, 2]
   -h, -?, --help          : This help message
   -m                      : Generate figures
   -p [plate_filename]     : Plate filepath
   --python [file path]    : Use specified Python executable [Default:
/usr/local/bin/python3]
   -s [sample_filename]    : Sample name and replicate file [used with -t 3 or 4]
   -t [input_file_type]    : Format number of PM plate format. Range: [1, 6] [Default: 1]
   -v                      : Verbose output
```

## *Cite*
*Cuevas, D. A., & Edwards, R. A. (2017).* ***PMAnalyzer: a new web interface for bacterial growth curve analysis.*** *Bioinformatics, 33(12), 1905–1906. (doi: [10.1093/bioinformatics/btx084](http://dx.doi.org/10.1093/bioinformatics/btx084))*
