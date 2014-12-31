#!/usr/bin/perl
# pmParser.pl
# Parsing script for GT Analyst multi-plate
# spectrophotometer
#
# Author: Daniel A Cuevas
# Created on 22 Nov. 2013
# Updated on 31 Dec. 2014

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
                            conc   => $conc
        };
    }
    return $plate;
}


###
# Input is the file path to the PM analyst file. Data is stored in the
# data hash-reference object
###
sub readData {
    my ($file, $data, $reps, $time) = @_;
    my ($name, $rep);
    my $fname = fileparse($file, qr/\.[^.]*/); # Capture file extension

    # Grab name (and replicate if applicable)
    if( $reps ) {
        ($name, $rep) = $fname =~ /^([A-Za-z0-9-.]+)_([A-Za-z0-9]+)/;
    }
    # Set replicate to 1 if not applicable
    else {
        ($name) = $fname =~ /^([A-Za-z0-9-._]+)/;
        $rep = 1;
    }

    # Check that the name and replicate variables are defined
    ! defined $name && ! defined $rep ? &error("Could not extract name (and replicate) from $file") :
                        ! defined $name ? &error("Could not extract name from $file") :
                        ! defined $rep ? &error("Could not extract replicate from $file") : 1;

    # Begin data extraction
    open(FILE,"<$file") or die "Couldn't open $file for reading\n";

    # Time difference is calculated to record relative time
    # instead of current time
    my ($prevTime, $currTime, $bckgrnd);
    my @timePoints;
    while( <FILE> ) {
        chomp;
        # Check for background
        if( /BACKGROUND/ ) {
            $bckgrnd = 1;
        }
        # Find time point
        elsif( /(\d.*?\s.{11})\s+/ ) {
            # No longer in background
            $bckgrnd = 0;
            my $timeRead = str2time($1);
            if( ! defined $prevTime ) {
                push(@timePoints, "0.0");
                $currTime = "0.0";
            }
            else {
                my $diffTime = ($timeRead - $prevTime) / 3600;
                $currTime = sprintf("%.1f", $timePoints[-1] + $diffTime);
                push(@timePoints, $currTime);
            }
            # Set new previous time
            $prevTime = $timeRead;
        }
        elsif( ! $bckgrnd && /(\w\d+)\s+([0-9.]+)/ ) {
            # Well = $1
            # Data = $2
            #print STDERR "Well $1, Data $2\n";
            #<>;
            $data->{$name}->{$rep}->{$1}->{$currTime} = $2;
        }
    }
    @$time = @timePoints if @timePoints > @$time;
    close(FILE);
}


###
# Calculate the median of the replicates and return the new data
# hash. This hash will not have a replicates key.
###
sub calcMedian {
    use Statistics::Basic qw(median);
    my ($data) = @_;
    my $newdata = {};
    foreach my $c ( keys %$data ) {
        my $repCounter = 0;
        foreach my $r ( keys %{$data->{$c}} ) {
            ++$repCounter;
            foreach my $w ( keys %{$data->{$c}->{$r}} ) {
                my @ods = ();
                foreach my $t( keys %{$data->{$c}->{$r}->{$w}} ) {
                    my $val = ($data->{$c}->{$r}->{$w}->{$t}) ?
                                $data->{$c}->{$r}->{$w}->{$t} :
                                0;
                    if( $repCounter == 1 ) {
                        $newdata->{$c}->{$w}->{$t} = ();
                    }
                    push(@{$newdata->{$c}->{$w}->{$t}}, $val);
                }
            }
        }
    }
    foreach my $c ( keys %$newdata ) {
        foreach my $w ( keys %{$newdata->{$c}} ) {
            foreach my $t ( keys %{$newdata->{$c}->{$w}} ) {
                #print $t."\n";
                my @ods = @{$newdata->{$c}->{$w}->{$t}};
                my $numReps = @ods;
                my $medvector = median()->set_size($numReps);
                foreach( @ods ) {
                    $medvector->insert($_);
                }
                my $median = $medvector->query;
                $newdata->{$c}->{$w}->{$t} = $median;
            }
        }
    }
    return $newdata;
}


###
# Calculate the mean of the replicates and return the new data
# hash. This hash will not have a replicates key.
###
sub calcMean {
    use Statistics::Basic qw(mean);
    my ($data) = @_;
    my $newdata = {};
    foreach my $c ( keys %$data ) {
        my $repCounter = 0;
        foreach my $r ( keys %{$data->{$c}} ) {
            ++$repCounter;
            foreach my $w ( keys %{$data->{$c}->{$r}} ) {
                my @ods = ();
                foreach my $t( keys %{$data->{$c}->{$r}->{$w}} ) {
                    my $val = ($data->{$c}->{$r}->{$w}->{$t}) ?
                                $data->{$c}->{$r}->{$w}->{$t} :
                                0;
                    if( $repCounter == 1 ) {
                        $newdata->{$c}->{$w}->{$t} = ();
                    }
                    push(@{$newdata->{$c}->{$w}->{$t}}, $val);
                }
            }
        }
    }
    foreach my $c ( keys %$newdata ) {
        foreach my $w ( keys %{$newdata->{$c}} ) {
            foreach my $t ( keys %{$newdata->{$c}->{$w}} ) {
                #print $t."\n";
                my @ods = @{$newdata->{$c}->{$w}->{$t}};
                my $numReps = @ods;
                my $meanvector = mean()->set_size($numReps);
                foreach( @ods ) {
                    $meanvector->insert($_);
                }
                my $mean = $meanvector->query;
                $newdata->{$c}->{$w}->{$t} = $mean;
            }
        }
    }
    return $newdata;
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
    my ($data, $plate, $pn, $reps, $time) = @_;
    foreach my $c ( keys %$data ) {
        foreach my $r ( map { $_->[0] } sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_,/([A-Za-z]+)/,/(\d+)/] } keys %{$data->{$c}} ) {
            # If the replicates flag is set then we are in the wells already
            if( $reps ) {
                my $w = $r;
                my @ods = ();
                foreach my $t( @$time ) {
                    my $val = ($data->{$c}->{$w}->{$t}) ?
                                $data->{$c}->{$w}->{$t} :
                                0;
                    push(@ods, $val);
                }
                my $ms = $plate ? $plate->{$w}->{main_source} : 0;
                my $s = $plate ? $plate->{$w}->{substrate} : 0;
                $plate ?
                    print join("\t", ($c, $ms, $s, $pn, $w, @ods))."\n" :
                    print join("\t", ($c, $w, @ods))."\n"
                ;
            }
            else {
                foreach my $w ( map { $_->[0] } sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_,/([A-Za-z]+)/,/(\d+)/] } keys %{$data->{$c}->{$r}} ) {
                    my @ods = ();
                    foreach my $t( @$time ) {
                        my $val = ($data->{$c}->{$r}->{$w}->{$t}) ?
                                    $data->{$c}->{$r}->{$w}->{$t} :
                                    0;
                        push(@ods, $val);
                    }

                    my $ms = $plate ? $plate->{$w}->{main_source} : 0;
                    my $s = $plate ? $plate->{$w}->{substrate} : 0;
                    $plate ?
                        print join("\t", ($c, $ms, $s, $pn, $w, @ods))."\n" :
                        print join("\t", ($c, $w, @ods))."\n"
                    ;
                }
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
    &usage("Directory empty given") if $#filepaths == 0;
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
        &usage("Plate file does not exist");
    }
}


######################################################
## Data extraction process
######################################################
my $data = {};
my @time;

# Read in data for each file
foreach my $f ( @filepaths ) {
    &readData($f, $data, $opts->{reps}, \@time);
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
foreach my $timeIter( @time ) {
    print "\t".sprintf("%.1f", $timeIter);
}
print "\n"; # End of header line

# Print out data
&printData($data, $plate, $plateName, $opts->{reps}, \@time);
