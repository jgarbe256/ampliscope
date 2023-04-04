#!/usr/bin/perl -w

#####################################
# dedup.pl
# John Garbe
# June 2023
#####################################

=head1 DESCRIPTION

dedup.pl - deduplicate reads in a fastq file

=head1 SYNOPSIS

dedup.pl [--minreads 1] R1.fastq R2.fastq

=cut

# read in a fastq file, dedup, and spit out a fasta file

use Getopt::Long;
use Pod::Usage;

$args{minreads} = 1;
GetOptions("help" => \$args{help},
	   "verbose" => \$args{verbose},
	   "minreads=i" => \$args{minreads},
    ) or pod2usage;
pod2usage(-verbose => 99, -sections => [qw(DESCRIPTION|SYNOPSIS|OPTIONS)]) if ($args{help});

$ifile = shift @ARGV;
open IFILE, $ifile or die "cannot open input file $ifile: $!\n";

# read in sequences from fastq file
while ($line = <IFILE>) {
    $seq = <IFILE>;
    chomp $seq;
    $sequences{$seq}++;
    $line = <IFILE>;
    $line = <IFILE>;
}

$count = 0;
$total = 0;
foreach $seq (keys %sequences) {
    $total++;
    next if ($sequences{$seq} < $args{minreads});
    $count++;
    print ">$count:$sequences{$seq}\n$seq\n";
    
}
print STDERR "$count sequences with at least $args{minreads} reads printed\n";
$skipped = $total - $count;
print STDERR "$skipped sequences with less than $args{minreads} reads discarded\n";

$ofile = "$ifile.dedup.counts";
open OFILE, ">$ofile" or die "cannot open output file $ofile: $!\n";
foreach $seq (sort { $sequences{$a} <=> $sequences{$b} } keys %sequences ) {
    $counts{$sequences{$seq}}++;
}
foreach $i (sort {$a <=> $b } keys %counts) {
    print OFILE "$i\t$counts{$i}\n";
}

