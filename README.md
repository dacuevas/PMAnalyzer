##*What is PMAnalyzer?*
PMAnalyzer is a command-line software pipeline that performs bacterial growth curve analysis for **M**ulti-phenotype **A**ssay **P**lates (MAPs).
 MAPs are 96-well micro-titer plates used to measure biomass formation through a spectrophotometer. Through the use of Python's SciPy and NumPy packages,
 PMAnalyzer performs curve-fitting in order to assess lag phase, growth rate, biomass yield, and growth level on each growth curve.

##*Dependencies*
####Languages
1. Bash
2. Perl
3. Python 2.7.X+ or 3.X.X

####Modules/Packages
1. SciPy
2. NumPy
3. Matplotlib (for graphs and figures)

##*Input files*
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

##*Output files*
(*All output text files are tab-delimited*)
- Various growth curves text files (raw, mean, median, logistic)
- Growth curve parameter text files (raw, logistic)
- Various growth curve figures
- Growth level heat-map

####Sample Images
![heatmap](https://github.com/dacuevas/PMAnalyzer/blob/master/sample/sample_results/growthlevels.png "Growth Level Heatmap")

![growthcurves](https://github.com/dacuevas/PMAnalyzer/blob/develop/sample/sample_results/avg_R.S.3.png "C. sedlakii Growth Curves")

##*Usage*
```
        ___          _               _
       / _ \/\/\    /_\  _ __   __ _| |_   _ _______ _ __
      / /_)/    \  //_\\| '_ \ / _` | | | | |_  / _ \ '__|
     / ___/ /\/\ \/  _  \ | | | (_| | | |_| |/ /  __/ |
     \/   \/    \/\_/ \_/_| |_|\__,_|_|\__, /___\___|_|
                                        |___/

PMAnalyzer version 1.2

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
   --debug                 : Print out debug messages
   -e [encoding]           : Check for file encoding of original data [Default:
UTF-16]
   -f                      : Flag if filter should be applied to growth curves
during PManalysis
   -g                      : Flag to use new growth level calculation
   -h, -?, --help          : This help message
   -m                      : Generate figures
   -n [number of plates]   : Number of plates [used with some parsers]
   -p [plate_filename]     : Plate filepath
   --python [file path]    : Use specified Python executable [Default:
/usr/local/bin/python3.4]
   -s [sample_filename]    : Sample name and replicate file [used with some
parsers]
   -t [input_file_type]    : Format number of PM plate format. Range: [1, 2]
   -v                      : Verbose output
```

##*Cite*
*Cuevas DA, Garza D, Sanchez SE et al.* ***[Elucidating genomic gaps using phenotypic profiles](http://f1000research.com/articles/3-210/)***
 *[v1; ref status: approved with reservations 1, http://f1000r.es/488] F1000Research 2014, 3:210
 (doi: [10.12688/f1000research.5140.1](http://dx.doi.org/10.12688/f1000research.5140.1))*
