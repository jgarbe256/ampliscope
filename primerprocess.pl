#!/usr/bin/perl -w

#######################################################
# primerprocess.pl
# John Garbe
# February 2013
#
# Read in a deduplicated fasta file with sequence counts in the read name
# Read in a file with a list of primers and sample names
# Count the frequency of each sequence seen, for each sample
#
# Primer.txt format:
# sample_name forward_primer reverse_primer expectedlength
# note: reverse_primer should be the reverse complement of the actual reverse primer sequence
#
# TODO:
#
#######################################################

=head1 DESCRIPTION

primerprocess.pl

=head1 SYNOPSIS

primerprocess.pl primers.txt sample.fasta

Options:

 --help : Print usage instructions and exit
 --verbose : Print more information while running
 --threads 10 : Number of threads/cores to use
 --margin 10 : Keep sequences within +/- margin of the expected lenth
 --minreads 1 : Keep sequences with at least minreads reads
 --revcomp : Reverse complement the second column of primer sequences
 --primtlimit 1000 : # how many unique inserts to print out 
 --dimerlength 10 : # sequences shorter than the primer lengths plus dimerlength is considered a primer dimer
=cut


use Getopt::Long;
use Pod::Usage;

$args{minreads} = 1; # used 4 for horizon data
$args{printlimit} = 2000; 
$args{dimerlength} = 10; 
$args{threads} = 8;
GetOptions("help" => \$args{help},
	   "verbose" => \$args{verbose},
	   "minreads=i" => \$args{minreads},
	   "threads=i" => \$args{threads},
	   "revcomp" => \$args{revcomp},
	   "margin=i" => \$args{margin},
	   "printlimit=i" => \$args{printlimit},
	   "dimerlength=i" => \$args{dimerlength},
    ) or pod2usage;
pod2usage(-verbose => 99, -sections => [qw(DESCRIPTION|SYNOPSIS|OPTIONS)]) if ($args{help});
pod2usage unless ($#ARGV == 1);

$args{primerfile} = shift @ARGV;
$args{fastafile} = shift @ARGV;

# read in primer file
open PFILE, $args{primerfile} or die "Cannot open primer file $args{primerfile}\n";

print STDERR "Reading primer file...\n";
$primerlinecount = 0;
while (<PFILE>) {
    $primerlinecount++;
    next if ($_ =~ /^#/);
    chomp;
    next if ($_ eq "");
    @line = split;
    $columns = $#line+1 unless (defined($columns));

    die "Error: Primer file line $primerlinecount has " . $#line+1 . " columns, expecting $columns columns: $_\n" if ($#line+1 != $columns);

    ($samplename, $f1, $r1, $length) = @line;
    $r1 = &reverse_complement($r1) if ($args{revcomp});
    $barcodes{$samplename}{f1} = $f1 // "";
    $barcodes{$samplename}{r1} = $r1 // "";
    $barcodes{$samplename}{length} = $length - length($f1) - length($r1);
    if ($columns == 6) {
	$barcodes{$samplename}{start} = $line[4];
	$barcodes{$samplename}{end} = $line[5];
    }
}
$cutout = 0;
$cutout = 1 if ($columns == 6);

# read in fasta file
open IFILE, $args{fastafile} or die "cannot open fasta file $args{fastafile}: $!\n";
$fastaseqcount = `grep "^>" $args{fastafile} | wc -l`;
chomp $fastaseqcount;

print STDERR "Reading fasta file with $fastaseqcount reads...\n";
$totalseqcount = 0;
$totalreadcount = 0;
$matchseqcount = 0;
$matchreadcount = 0;
$primerdimerreadcount = 0;
$primerdimerseqcount = 0;
$progressincrement = 10;
$progresspct = $progressincrement;
$filterreadcount = 0;
$filterseqcount = 0;
$badlengthreadcount = 0;
$badlengthseqcount = 0;
while ($id = <IFILE>) {
    if ($totalseqcount / $fastaseqcount * 100 > $progresspct) {
	print STDERR "$totalseqcount sequences processed ($progresspct%)...\n";
	$progresspct += $progressincrement;
    }
    $seq = <IFILE>;
    chomp $id;
    my ($junk, $seqcount) = split /:/, $id;
    $totalseqcount++;
    $totalreadcount += $seqcount;
    if ($seqcount < $args{minreads}) {
	$filterreadcount += $seqcount;
	$filterseqcount++;
	next;
    }
    chomp $seq;
    $qes = reverse_complement($seq);
    $match = 0;
    $matchedseq = "";
    $matchedinsert = "";
    for $barcode ( keys %barcodes ) {
	$f1 = $barcodes{$barcode}{f1};
	$r1 = $barcodes{$barcode}{r1};
	if ($seq =~ /^$f1([agctnAGCTN]+)$r1/i) {
	    $matchedseq = $seq;
	    $matchedinsert = $1;
	} elsif ($qes =~ /^$f1([agctnAGCTN]+)$r1/i) { # handle reverse complements
	    $matchedseq = $qes;
	    $matchedinsert = $1;
	}
	if ($matchedseq) {
	    $totalreads{$barcode} += $seqcount;
	    if (length($matchedinsert) < $args{dimerlength}) {
		# primerdimer
		$primerdimerseqcount{$barcode}++;
		$primerdimerreadcount{$barcode} += $seqcount;
		$primerdimerseqcount++;
		$primerdimerreadcount += $seqcount;
		next;
	    }
	    if ((length($matchedinsert) < ($barcodes{$barcode}{length} - $args{margin})) or (length($matchedinsert) > ($barcodes{$barcode}{length} + $args{margin}))) {
		# too short or too long;
		$badlengthseqcount{$barcode}++;
		$badlengthreadcount{$barcode} += $seqcount;
		$badlengthseqcount++;
		$badlengthreadcount += $seqcount;
		next;
	    }
	    if ($cutout) { # cut out a region of interest in the amplicon
		$front = $barcodes{$barcode}{start};
		$back = $barcodes{$barcode}{end};
		if (length($matchedseq) < $back) {
#		print "too short: " . length($qes) . " $front $back\n";
		    next;
		}
		$end = $back - $front;
		$insert = substr($matchedseq,$front,$end);
	    } else {
		$insert = $seq;
	    }
	    if ($inserts{$barcode}{$insert}) { # DEBUG
#		print "insert seq already seen: $insert\n";
	    }
	    $inserts{$barcode}{$insert} += $seqcount; # debug, was = 
	    $matchseqcount{$barcode}++;
	    $matchreadcount{$barcode} += $seqcount;
	    $matchseqcount++;
	    $matchreadcount += $seqcount;
	    $match = 1;
	    last; # barcode found, move on to next sequence

	}
    }
    if ($match == 0) { # keep track of what we aren't matching
	$nobarcodeseqcount++;
	$nobarcodereadcount += $seqcount;
	$bad = substr($seq, 0, 20) . "-" . substr($seq, -20);
	if (defined($badbarcodes{$bad})) {
	    $badbarcodes{$bad} += $seqcount;
	} else {
	    $dab = &reverse_complement($bad);
	    if (defined($badbarcodes{$dab})) {
		$badbarcodes{$dab} += $seqcount;
	    } else {
		$badbarcodes{$bad} = $seqcount;
	    }
	}
    }
}

print STDERR "$totalseqcount seqs, $totalreadcount reads\n";
print STDERR "$nobarcodeseqcount nobarcodeseqs, $nobarcodereadcount nobarcodereads\n";
print STDERR "$badlengthseqcount badlengthseqs, $badlengthreadcount badlengthreads\n";
print STDERR "$filterseqcount filterseqs, $filterreadcount filterreads\n";
print STDERR "$primerdimerseqcount primerdimerseqs, $primerdimerreadcount primerdimerreads\n";
print STDERR "$matchseqcount matchseqs, $matchreadcount matchreads\n";

#$percentmatch = int( ($matchcount / $count) * 1000) / 10;
#print STDERR "\n$count reads processed\n";
#print STDERR "$matchcount reads matched to barcode ($percentmatch%)\n";

if ($matchreadcount / $totalreadcount < .05) {
    print STDERR "Warning: few reads matched primer sequences. Check primer sequences and expected lengths, try using the --revcomp option\n";
}
die "Error: no sequences remaining after deduplication\n" if ($matchreadcount == 0);

# print out demultiplexing stats
$ofile = "demultiplex-stats.txt";
open OFILE, ">$ofile" or die $!;
print OFILE "Amplicon\tTotalReads\tPassedSequences\tPassedReads\tDimerReads\tBadLengthReads\n";
for $barcode ( sort keys %barcodes ) {
    $matchseqcount = $matchseqcount{$barcode} // 0; #insertcount = 0;
    $matchreadcount = $matchreadcount{$barcode} // 0; #seqcount = 0;
#    for $insert ( sort keys %{ $inserts{$barcode} } ) {
#	$insertcount++;
#	$seqcount += $inserts{$barcode}{$insert};
#    }
    $badlength = $badlengthreadcount{$barcode} // 0;
    $primerdimer = $primerdimerreadcount{$barcode} // 0;
    $totalreads = $totalreads{$barcode} // 0;
    print OFILE "$barcode\t$totalreads\t$matchseqcount\t$matchreadcount\t$primerdimer\t$badlength\n";
}
if ($args{verbose}) {
    $text = `cat $ofile`;
    print STDERR $text;
}
#    exit; # uncomment to end here

# print out the unique sequences seen to fasta file
for $barcode ( sort keys %barcodes ) {
    $insertcount = 0;
    $addnewline = 0;
    $othercount = 0;
    $ofile = $barcode . ".demux.fa";
    open OFILE, ">$ofile" or die "Cannot open output file $ofile\n";
#    print STDERR "amplicon\tinsert length\tinsert count\tinsert sequence\n";
    for $insert ( sort { $inserts{$barcode}{$b} <=> $inserts{$barcode}{$a} } keys %{ $inserts{$barcode} } ) {
	$insertcount++;
	$addnewline = 1;
#	$insertlength= length($insert);
	if ($insertcount > $args{printlimit}) {
	    $othercount += $inserts{$barcode}{$insert};
	} else {
#	    print STDERR "$barcode\t$insertlength\t$inserts{$barcode}{$insert}\t$insert\n";
	    print OFILE ">$barcode:$inserts{$barcode}{$insert}\n$insert\n";
	}
    }
#    print STDERR "$barcode\tNA\t$othercount\tother inserts\n" unless ($othercount == 0);
#    print STDERR "\n" if ($addnewline == 1);
    close OFILE;
}

# print out barcodes that weren't matched
$limit = 40;
$ofile = "unmatched-barcodes.txt";
open OFILE, ">$ofile" or die "$!\n";
print OFILE "\n\n### Unmatched barcodes (top $limit) ###\n";
print OFILE "barcode\tReads\n";
$count = 0;
for $barcode ( sort { $badbarcodes{$b} <=> $badbarcodes{$a} } keys %badbarcodes ) {
    print OFILE "$barcode\t$badbarcodes{$barcode}\n";
    $count++;
    last if ($count > $limit); # limit how many unmatched are printed
}
close OFILE;
$text = `cat $ofile`;
print STDERR $text if ($args{verbose});

########################### helper subs #################################

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
#    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}
