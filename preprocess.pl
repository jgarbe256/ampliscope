#!/usr/bin/perl -w

##################################################
# preprocess.pl
# John Garbe
# August 2016
#
# Pre-process amplicon data:
#  Merge paired-reads
#  Deduplicate (and count)
#
###################################################

use File::Basename;
use Cwd 'abs_path';
use Getopt::Std;

$usage = "USAGE: preprocess.pl R1.fastq [ R2.fastq ] [minlength] [maxlength]\n";
die $usage unless getopts('h');
die $usage unless (($#ARGV >= 0) and ($#ARGV <= 3));

my ($r1, $r2, $min, $max) = @ARGV;

our ($opt_h);
die $usage if ($opt_h);

$THREADS = 40;

my ($samplename) = fileparse($r1);
print "$samplename\n";
$samplename =~ /(.*)_R1/;
$name = $1 // "";
if (! $name) {
    $samplename =~ /(.*)\.fast[qa]/;
    $name = $1 // "";
}
$samplename = $name || "sample";
print STDERR "processing sample $samplename\n";

die "Cannot find file r1 file $r1\n" unless (-e $r1);
$r1 = abs_path($r1);
if ($r2) {
    die "Cannot find file r2 file $r2\n" unless (-e $r2);
    $r2 = abs_path($r2);
}

### merge reads with pear ###
print STDERR "Merging reads with pear...\n";
if ($r2) {
    $length = "";
    if ($min and $max) {
	$length = "-n $min -m $max";
    } else {
	$length = "-n 30"; # this will help stitch primer dimers
    }
    `pear -f $r1 -r $r2 -e -o $samplename $length -j $THREADS > $samplename.pear.log`; # don't calculate empirical frequencies: calculation takes a long-ass time
    $fastqfile = "$samplename.assembled.fastq";
} else {
    print STDERR "single-end reads, skipping merging\n";
    $fastqfile = $r1;
}

# print out # and % of reads merged
$result = `grep "^Assembled reads" $samplename.pear.log`;
print $result;

### Dedup with dedup.pl ###
print STDERR "Deduplicating reads with dedup.pl...\n";
$deduped = "$samplename.deduped.fasta";
`dedup.pl $fastqfile > $deduped`;

exit;
