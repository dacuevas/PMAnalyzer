#!/usr/bin/perl
# pmParser_v2.pl
# Parsing script for Molecular Devices multi-plate
# spectrophotometer.
#
# Author: Daniel A Cuevas
# Created on 18 Mar. 2015
# Updated on 02 Apr. 2015
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file of continous save feature
# === Each plate recording is appended each time a read occurs
#     Example:
#              ##BLOCKS= 1
#              Plate:   Plate1
#                       Temperature    A1      A2
#                       25.0          0.123    0.135
#
#               ~End
#               Original filename: ... Date Last Saved: 2/27/2015 2:00:00PM
#              ##BLOCKS= 1
#              Plate:   Plate1
#                       Temperature    A1      A2
#                       25.0          0.123    0.135
#
#               ~End
#               Original filename: ... Date Last Saved: 2/27/2015 2:00:00PM
#

use warnings;
use strict;
use Date::Parse;
use File::Basename;
use Getopt::Long;


################################################
## Subroutines
################################################

###
# Printout when command line arguments are incorrect
###
sub usage {
    my $msg = shift;
    my $filepath = basename($0);
    print "ERROR: $msg\n\n" if $msg;
    print << "END"
USAGE
    perl $filepath [Options] <PM_file.txt OR PM_directory>

REQUIRED
    -f, --files "Data file list"    : Comma-separated list of data files
               OR
    -d, --dir "Data directory"      : Directory containing data files
                                       (Given preference over -f)
    -n, --nplates "No. plates"      : Number of plates
    -s, --samples "samples.txt"     : Tab-delimited file of sampe names

OPTIONS
    -p, --plate "platePath.txt"     : Plate filepath
    -?, -h, --help                  : This usage message

NOTES
    FILE NAME FORMATS:
        When using replicates be sure to have the file name in the format,

        <Sample Name>_<Replicate Letter>_<other text>.txt
        Regex used:    <Sample Name> -- [A-Za-z0-9-.]+
                       <Replicate Letter> -- [A-Za-z0-9]+

        When replicates are not indicated,
        <Sample Name>_<other text>.txt
        Regex used:    <Sample Name> -- [A-Za-z0-9-._]+

END
;
    exit(1);
}


###
# Error message output and exit
###
sub error {
    my $msg = shift;
    print STDERR "ERROR: $msg\n";
    exit(2);
}


###
# Input is array of file paths. Files are checked for existence
###
sub checkFiles {
    my $files = shift;
    foreach( @$files ) {
        return 0 if ! -e $_;
    }
    return 1;
}


###
# Input is the file path to the plate file. Plate information is stored
# (and returned) as a hash-reference object
###
sub readPlate {
    my $file = shift;
    open(FILE,"<$file") or die "Couldn't open $file for reading.\n";
    my $plate = {};
    while( <FILE> ) {
        chomp;
        my ($well, $mainSource, $substrate, $conc) = split /\t/;
        $plate->{$well} = { main_source => $mainSource,
                            substrate  => $substrate,
                            conc   => $conc };
    }
    return $plate;
}


###
# Input is the file path to the sample name file. Sample name and replicate
# information is stored (and returned) as a hash-reference object
###
sub readSamples {
    my $file = shift;
    open(FILE,"<$file") or die "Couldn't open $file for reading.\n";
    my $samples->[0] = 0;
    my $idx = 0;
    while( <FILE> ) {
        chomp;
        my ($name, $rep) = split /\t/;
        $samples->[$idx++] = { name => $name,
                               rep  => $rep };
    }
    return $samples;
}


###
# Input is the file path to the PM analyst file. Data is stored in the
# data hash-reference object
###
sub readData {
    my ($file, $data, $time, $t0, $wells, $nPlates, $sNames) = @_;

    # Begin data extraction
    open(FILE,"<$file") or die "Couldn't open $file for reading\n";

    # Time difference is calculated to record relative time
    # instead of current time
    my ($currTime, $dataFound, @dataPoints, @timePoints, $name, $rep);
    my $plateIdx = 0;
    while( <FILE> ) {
        chomp;
        next if /^$/;  # Skip blank lines
        # Check if we find data line
        my @l = split(/\t/);
        if( $#l > 90 && $_ !~ /[A-z]/ ) {
            @dataPoints = @l;
            # At the data line
            # Remove blank items
            my $shifted = shift @dataPoints;
            # Remove temperature
            $shifted = shift @dataPoints;
            # Get sample name and replicate
            $name = $sNames->[$plateIdx]{name};
            $rep = $sNames->[$plateIdx]{rep};
        }
        elsif( /(\d{1,2}\/\d{1,2}\/\d{4} \d{1,2}:\d{2}:\d{2} \w{2})/ ) {
            # Get time from line
            my $recTime = str2time($1);
            my $currTime;
            # Store time difference
            if( !defined $time->{$name}->{$rep} ) {
                $time->{$name}->{$rep} = [];
                $t0->{$name}->{$rep} = $recTime;
                $currTime = "0.0";
            }
            else {
                $currTime = ($recTime - $t0->{$name}->{$rep}) / 3600;
                $currTime = sprintf("%.1f", $currTime);
            }
            push(@{$time->{$name}->{$rep}}, $currTime);

            foreach my $idx ( 0..$#dataPoints ) {
                # Should be 96
                my $w = $wells->[$idx];
                my $d = $dataPoints[$idx];
                if( !defined $data->{$name}->{$rep}->{$w} ) {
                    $data->{$name}->{$rep}->{$w} = [$d];
                }
                else {
                    push(@{$data->{$name}->{$rep}->{$w}}, $d);
                }
            }
            # Reset for next plate info
            $plateIdx = ($plateIdx + 1) % $nPlates;
        }
    }
    close(FILE);
}


###
# Input is the data hash-reference object. Information is printed in the order:
# (All output is tab-delimited)
# Row 1 (header) -- sample (main_source substrate well# plate_name) [time]
# Row 2-n -- sample (main_source substrate well# plate_name) [data]
# Note: When replicates are kept separate, the sample is separated into two
# columns for "sample" and "rep"
###
sub printData {
    my ($data, $plate, $pn, $time) = @_;
    my $timeLen = scalar @$time;
    foreach my $c ( keys %$data ) {
        foreach my $r ( map { $_->[0] } sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_,/([A-Za-z]+)/,/(\d+)/] } keys %{$data->{$c}} ) {
            foreach my $w ( map { $_->[0] } sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_,/([A-Za-z]+)/,/(\d+)/] } keys %{$data->{$c}->{$r}} ) {
                my @ods = ();
                foreach my $tIdx ( 0..$timeLen-1 ) {
                    my $val = ($data->{$c}->{$r}->{$w}->[$tIdx]) ?
                                $data->{$c}->{$r}->{$w}->[$tIdx] :
                                0;
                    push(@ods, $val);
                }

                my $cr = "${c}_${r}";
                my $ms = $plate ? $plate->{$w}->{main_source} : 0;
                my $s = $plate ? $plate->{$w}->{substrate} : 0;
                $plate ?
                    print join("\t", ($cr, $ms, $s, $pn, $w, @ods))."\n" :
                    print join("\t", ($cr, $w, @ods))."\n"
                ;
            }
        }
    }
}


###
# Custom sorting method for files
###
sub fileSorter {
    my ($a_name, $a_rep) = $a =~ /([A-z0-9-]+)[_\s](\d+).txt$/;
    my ($b_name, $b_rep) = $b =~ /([A-z0-9-]+)[_\s](\d+).txt$/;
    return $a_name cmp $b_name || $a_rep <=> $b_rep;
}

###################################################
## ARGUMENT PARSING
###################################################

my $opts = {
                    files    =>   "",
                    dir      =>   "",
                    mean     =>   "",
                    nplates  =>   "",
                    plate    =>   "",
                    samples  =>   "",
                    help     =>   ""
};

GetOptions(
                    "f=s"          =>   \$opts->{files},
                    "files=s"      =>   \$opts->{files},
                    "d=s"          =>   \$opts->{dir},
                    "dir=s"        =>   \$opts->{dir},
                    "m"            =>   \$opts->{mean},
                    "mean"         =>   \$opts->{mean},
                    "n=s"          =>   \$opts->{nplates},
                    "nplates=s"    =>   \$opts->{nplates},
                    "p=s"          =>   \$opts->{plate},
                    "plate=s"      =>   \$opts->{plate},
                    "s=s"          =>   \$opts->{samples},
                    "samples=s"    =>   \$opts->{samples},
                    "h"            =>   \$opts->{help},
                    "help"         =>   \$opts->{help}
);

&usage("") if $opts->{help} || !$opts->{nplates} || !$opts->{samples} ||
              (!$opts->{files} && !$opts->{dir});


#####################################################
## Check if data files exist
#####################################################
my @filepaths;
if( $opts->{files} ) {
    @filepaths = split(",", $opts->{files});
}
# Check if directory exists and is not empty
else {
    &usage("Data directory does not exist") unless ! $opts->{dir} || -e  $opts->{dir};
    opendir(DIR, $opts->{dir}) or &usage("Couldn't open $opts->{dir} to read files");
    while( my $f = readdir(DIR) ) {
        # Skip directories beginning with "\."
        next if $f =~ m/^\./;
        push(@filepaths, "$opts->{dir}/$f");
    }
    # Check if directory was empty
    &usage("Empty directory given") if $#filepaths == -1;
}

# Check that files actually exist
&checkFiles(\@filepaths) || &usage("File(s) does not exist");

# Check if plate file exists
my $plate;
my $plateName;
if( $opts->{plate} ) {
    if( -e $opts->{plate} ) {
        $plate = &readPlate($opts->{plate});
        $plateName = fileparse($opts->{plate}, qr/\.[^.]*/); # Capture file extension
    }
    else {
        &usage("Plate file not found");
    }
}

# Check if samples file exists
my $samples;
if( -e $opts->{samples} ) {
    $samples = &readSamples($opts->{samples});
}
else {
    &usage("Samples file not found");
}


######################################################
## Data extraction process
######################################################
my $data = {};
my $time = {};
my $t0 = {};  # Starting time for each sample
my $wells = [];
foreach my $w1 ("A".."H") {
    foreach my $w2 (1..12) {
        push(@$wells, "$w1$w2");
    }
}

# Read in data for each file
my @sortFs =  sort { fileSorter } @filepaths;
foreach my $f ( @sortFs ) {
    &readData($f, $data, $time, $t0,
              $wells, $opts->{nplates}, $samples);
}

# If replicates exist, calculate median or mean
#$data = !$opts->{reps}  ? $data
#        : $opts->{mean} ? &calcMean($data)
#        :                 &calcMedian($data);

######################################################
## Data printout process
######################################################

# Print out headers first
# Depends on if a plate file was supplied or not
$plate ?
    print "sample\tmainsource\tsubstrate\tplate\twell" :
    print "sample\twell";

# Print out hours
# Just use first sample's time
my @samples = keys %$time;
my @sreps = keys %{$time->{$samples[0]}};
my $ptime = $time->{$samples[0]}->{$sreps[0]};
foreach my $timeIter( @$ptime ) {
    print "\t".sprintf("%.1f", $timeIter);
}
print "\n";

# Print out data
&printData($data, $plate, $plateName, $ptime);
