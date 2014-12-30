##*What is PMAnalyzer?*
PMAnalyzer is a command-line software pipeline that performs bacterial growth curve analysis for **M**ulti-phenotype **A**ssay **P**lates (MAPs).
 MAPs are 96-well micro-titer plates used to measure biomass formation through a spectrophotometer. Through the use of Python's SciPy and NumPy packages,
 PMAnalyzer performs curve-fitting in order to assess lag phase, growth rate, biomass yield, and growth level on each growth curve.

##*Dependencies*
####Languages
1. Bash
2. Perl
3. Python 3.X.X

####Modules/Packages
1. SciPy
2. NumPy
3. Matplotlib (for graphs and figures)

##*Input files*
The parsing scripts from PMAnalyzer are generated from the Molecular Devices Analyst GT multi-plate plate reader. Sample input files can be viewed in the
 `sample/data_csedlakii` directory.

##*Output files*
(*All output text files are tab-delimited*)
- Various growth curves text files (raw, mean, median, logistic)
- Growth curve parameter text files (raw, logistic)
- Various growth curve figures
- Growth level heat-map

##*Cite*
*Cuevas DA, Garza D, Sanchez SE et al.* ***[Elucidating genomic gaps using phenotypic profiles](http://f1000research.com/articles/3-210/)***
 *[v1; ref status: approved with reservations 1, http://f1000r.es/488] F1000Research 2014, 3:210
 (doi: [10.12688/f1000research.5140.1](http://dx.doi.org/10.12688/f1000research.5140.1))*
