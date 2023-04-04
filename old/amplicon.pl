#!/usr/bin/perl -w

##################################################
# ampliscope.pl
# John Garbe
# August 2016
#
# Process amplicon data:
#  Group reads by barcode
#  Plot most frequent reads per barcode
#  Generate plot showing indels and mismatches at each base of each barcode
#
###################################################

use File::Basename;
use Cwd 'abs_path';
use Getopt::Std;

$usage = "USAGE: amplicon.pl primers.txt sample.fastq\n";
die $usage unless getopts('hd');
die $usage unless ($#ARGV == 1);

my ($primerfile, $fastq) = @ARGV;

our ($opt_h, $opt_d);
die $usage if ($opt_h);
#pod2usage(q(-verbose) => 3) if ($opt_h);

$THREADS = 40;

my ($samplename) = fileparse($fastq);
print "$samplename\n";
$samplename =~ /(.*)_R1/;
$name = $1 // "";
if (! $name) {
    $samplename =~ /(.*)\.fast[qa]/;
    $name = $1 // "";
}
$samplename = $name || "sample";
print STDERR "processing sample $samplename\n";

die "Cannot find primer file $primerfile\n" unless (-e $primerfile);
$primerfile = abs_path($primerfile);
die "Cannot find file fastq file $fastq\n" unless (-e $fastq);
$fastq = abs_path($fastq);

`mkdir amplicon-$samplename`;
chdir "amplicon-$samplename";

### Segregate by primer ###
print STDERR "Primer processing...\n";
`primerprocess.pl $primerfile $fastq > $samplename.txt`;
#exit; # debug
### Generate gnuplots ###
print STDERR "Generating gnuplots with align-parse.pl...\n";
$commandfile = "align-parse-commands.txt";
open OFILE, ">$commandfile" or die "cannot open $commandfile: $!\n";
@fas = `ls *.fa`;
foreach $fa (@fas) {
    chomp $fa;
    print OFILE "align-parse.pl -v $fa\n";
}
close OFILE;
`cat $commandfile | parallel -j $THREADS`;
#`find *.fa -exec align-parse.pl {} \\;`;

### Generate html index ###
print STDERR "Generating html index...\n";
$ofile = "$samplename.html";
open OFILE, ">$ofile" or die "cannot open output file $ofile: $!\n";
@fas = `ls *.fa`;
foreach $fa (sort hdsort @fas) {
    chomp $fa;
    $name = $fa;
    $name =~ s/\.fa$//;
    $totalcount = `cut -f 2 $name.totalcount | head`;
    chomp $totalcount;
    $png = "$fa.png";
    print OFILE "<img src=\"png\\$png\"><br>\n";
    $html = "$fa.html";
    print OFILE "<a href=\"html\\$html\">$name</a> $totalcount reads<br>\n";
}

# log-scale version
$ofile = "$samplename-log.html";
open OFILE, ">$ofile" or die "cannot open output file $ofile: $!\n";
@fas = `ls *.fa`;
foreach $fa (sort hdsort @fas) {
    chomp $fa;
    $name = $fa;
    $name =~ s/\.fa$//;
    $totalcount = `cut -f 2 $name.totalcount | head`;
    chomp $totalcount;
    $png = "$fa-log.png";
    print OFILE "<img src=\"png\\$png\"><br>\n";
    $html = "$fa.html";
    print OFILE "<a href=\"html\\$html\">$name</a> $totalcount reads<br>\n";
}

# generate list of reads per amplicon
`cat *.totalcount > $samplename.totalcounts.txt`;

### Clean up ###
`mkdir png`;
`mv *.png png`;

`mkdir html`;
`mv *.fa.html html`;

`mkdir fa`;
`mv *.fa fa`;

`mkdir align`;
`mv *.aln align`;

print STDERR "Finished. Have a nice day!\n";


### sort function for horizon discovery
sub hdsort {
    my ($mya) = split /\./, $a;
    my ($myb) = split /\./, $b;
    return $mya cmp $myb;
}
