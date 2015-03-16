#!/usr/bin/perl
# well_vs_time_parser.pl
# Parsing script for Molecular Devices single-plate
# spectrophotometer.
#
# Author: Daniel A Cuevas
# Created on 16 Mar. 2015
# Updated on 16 Mar. 2015
#
#
###############################################################
#                       ASSUMPTIONS
###############################################################
# = Tab delimited file
# = First row is "Time" \t "Temperature" \t [well numbers]
# === Well numbers (and data) start on column 3
# === Well numbers are ordered as A1, A2,..., D1, D2,..., H12
# = First column are times
# === Time format is hh:mm:ss. Hour is optional
# === Regex: (\d{1,3})?:\d{2}:\d{2}
# = There is always extra text at the beginning and end of the document

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

OPTIONS
    -m, --mean                      : Calculate mean instead of median
    -p, --plate "platePath.txt"     : Plate filepath
    -r, --reps                      : Replicate flag
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
# Input is the file path to the PM analyst file. Data is stored in the
# data hash-reference object
###
sub readData {
    my ($file, $data, $reps, $time, $t0, $wells) = @_;
    my ($name, $rep);
    my $fname = fileparse($file, qr/\.[^.]*/); # Capture file extension

    # Grab name (and replicate if applicable)
    ($name, $rep) = $fname =~ /^([A-Za-z0-9-.]+)_([A-Za-z0-9-.]+)/;

    # Check that the name and replicate variables are defined
    ! defined $name && ! defined $rep ? &error("Could not extract name and replicate from $file") :
                       ! defined $name ? &error("Could not extract name from $file") :
                       ! defined $rep ? &error("Could not extract replicate from $file") : 1;

    # Begin data extraction
    open(FILE,"<$file") or die "Couldn't open $file for reading\n";

    # Time difference is calculated to record relative time
    # instead of current time
    my ($prevTime, $currTime, $dataFound);
    my @timePoints;
    while( <FILE> ) {
        chomp;
        if( !$dataFound && $_ !~ /^Time/ ) {
            # Skip lines until Time line is found
            next;
        }
        elsif( !$dataFound ) {
            # Time line is found here, so data comes after
            $dataFound = 1;
            next;
        }
        elsif( /^$/ || /End/ ) {
            last;
        }

        # At the data line now
        my @datList = split(/\t/);
        # Store time difference
        my $timeRead = shift @datList;
        # Make sure time is in YYYY-MM-DDThh:mm:ss format
        # Use the length to determine raw format
        my ($hour, $min, $sec);
        my $tLen = length($timeRead);
        if( $tLen == 4 || $tLen == 5 ) {
            # Format: m:ss or mm:ss
            ($min, $sec) = /(\d{1,2}):(\d{2})/;
            $hour = 0;
        }
        elsif( $tLen == 7 || $tLen == 8 ) {
            # Format: h:mm:ss or hh:mm:ss
            ($hour, $min, $sec) = /(\d{1,2}):(\d{2}):(\d{2})/;
        }
        # Convert time to seconds
        my $cTime = 3600 * $hour + 60 * $min + $sec;
        # Use 2001-01-01 at midnight as the "starting point"
        my $recTime = str2time("2001-01-01T00:00:00") + $cTime;

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

        shift @datList;  # Second value is temperature reading

        # Insert values into data hash
        my $numWells = scalar @$wells;  # Should be 96
        foreach my $idx ( 0..$numWells-1 ) {
            my $w = $wells->[$idx];
            my $d = $datList[$idx];
            if( !defined $data->{$name}->{$rep}->{$w} ) {
                $data->{$name}->{$rep}->{$w} = [$d];
            }
            else {
                push(@{$data->{$name}->{$rep}->{$w}}, $d);
            }
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

###################################################
## ARGUMENT PARSING
###################################################

my $opts = {
                    files    =>   "",
                    dir      =>   "",
                    mean     =>   "",
                    plate    =>   "",
                    reps     =>   "",
                    help     =>   ""
};

GetOptions(
                    "f=s"      =>   \$opts->{files},
                    "files=s"  =>   \$opts->{files},
                    "d=s"      =>   \$opts->{dir},
                    "dir=s"    =>   \$opts->{dir},
                    "m"        =>   \$opts->{mean},
                    "mean"     =>   \$opts->{mean},
                    "p=s"      =>   \$opts->{plate},
                    "plate=s"  =>   \$opts->{plate},
                    "r"        =>   \$opts->{reps},
                    "reps"     =>   \$opts->{reps},
                    "h"        =>   \$opts->{help},
                    "help"     =>   \$opts->{help}
);

&usage("") if $opts->{help} ||
            ( !$opts->{files} && !$opts->{dir} );


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
foreach my $f ( @filepaths ) {
    &readData($f, $data, $opts->{reps}, $time, $t0, $wells);
}

# If replicates exist, calculate median or mean
$data = !$opts->{reps}  ? $data
        : $opts->{mean} ? &calcMean($data)
        :                 &calcMedian($data);

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
