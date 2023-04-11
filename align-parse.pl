#!/usr/bin/perl -w

##############################################
# align-parse.pl
# John Garbe
# April 2013
#
# Given a multi-fasta file, align each sequence against the first, compile
# results (each sequence must be on one line)
#
# Given a pairwise blast output (outfmt 3), count the number of
# inserts, deletions, and mismatches at each base
#
#############################################

use Getopt::Std;

$usage = "USAGE: align-parse.pl sample.fa\n";
die $usage unless getopts('hvi:');
die $usage unless ($#ARGV == 0);
die $usage if ($opt_h);
$verbose = $opt_v // 0;

our ($opt_v, $opt_h, $opt_i);

die "USAGE: align-parse.pl [-i .95] inserts.fasta > counts.txt\n" unless ($#ARGV == 0);

$sensitivitylimit = 0.1; # smallest percent that gets drawn on log-scale plot

# set up reference sequence
$args{minidentity} = $opt_i // 0;
$ifile = $ARGV[0];
print STDERR "Processing $ifile\n" if ($verbose);
$ofile = $ifile;
$ofile =~ s/demux\.fa/good\.fa/;
$badfile = $ifile;
$badfile =~ s/demux\.fa/bad\.fa/;
$inserts = `wc -l < $ifile`;
chomp $inserts;
($inserts) = split ' ', $inserts; 
$inserts = $inserts / 2;
$afile = "$ifile.aln";
$reffa = "$ifile.ref";
$altfa = "$ifile.alt";
`head -n 2 $ifile > $reffa`;
$refid = `head -n 1 $reffa`;
chomp $refid;
$refid =~ s/^>//;
@line = split /:/, $refid;
$refcount = pop @line;
$reftrimid = $line[0];
$refseq = `tail -n 1 $reffa`;
chomp $refseq;
@refseq = split "", $refseq;
`cat $reffa > $ofile`;

# initialize counts to zero
#for $i (0..2) {
#    for $j (0..$#refseq) { # cludge, should figure out length of sequence
#	$count[$i][$j] = 0;
#    }
#}

# for each non-reference insert
$totalcount = $refcount;
$output = ($verbose) ? "" : "2> /dev/null";
for $i (1..$inserts-1) {
#for $i (1..1) { # DEBUG
    # align each insert to reference
    $line = $i * 2 + 2;
##    print STDERR "insert $i\tline $line\n";
    `head -n $line $ifile | tail -n 2 > $altfa`;
    $result = `needle $reffa $altfa -gapopen 10 -gapextend 0.5 $afile -aformat3 markx2 -endweight $output 2>&1`;
    if ($?) {
	print "needle error: $result\n";
    }

    # ignore sequences dissimilar to ref
    $identity = `grep Identity $afile`;
    $identity =~ /\((\d+\.\d+)\%\)/;
    if ($1 < $args{minidentity}) {
	`cat $altfa >> $badfile`;
	next 
    }
    `cat $altfa >> $ofile`;

    # grab the count from the ID line
    $header = `head -n 1 $altfa`;
    chomp $header;
    @line = split /:/, $header;
    $altcount = pop @line;
    $totalcount += $altcount;
##    print "count: $altcount\n";

    # pull the reference (subject) and alternative (query) sequence
    # out of the messy needle output
    $ref = "";
    $alt = "";
    open AFILE, "$afile" or die "cannot open blast output $afile\n";
    $needlerefid = substr($refcount,0,6); # this works when the delimiter is :
#    $needlerefid = substr($refid,0,6);
    while (<AFILE>) {
	next unless (/^ *$needlerefid /);
	@line = split ' ';
	$ref .= $line[1];
	$line = <AFILE>;
	chomp $line;
	@line = split ' ', $line;
	$alt .= $line[1];
    }

    print "reference: $ref\n";
    print "query    : $alt\n";

    @ref = split "", $ref;
    @alt = split "", $alt;

    # count each insertion, deletion, and mismatch
    $offset = 0; # inserts mess up the base coordinates, fix with an offset
#    for $j (0..$#ref) {
    $j = 0;
    while ($j <= $#ref) {
#	print "J: $j, ref: $ref[$j], alt: $alt[$j]\n";
	if ($ref[$j] eq "-") { # insertion
	    # consume insertion
	    $insertion = "$alt[$j]";
	    while (($j < $#ref) and ($ref[$j+1] eq "-")) {
		$j++;
		$offset++;
		$insertion .= "$alt[$j]";
	    }
#	    print "J: $j, offset: $offset, insert: $insertion\n";
#	    print "insertion: $insertion ($j $offset)\n";
	    $count[$j-$offset]{in}{$insertion} += $altcount;
	    $count[$j-$offset]{intotal} += $altcount;
	    $offset++;
	} elsif ($alt[$j] eq "-") { # deletion
	    # consume deletion
	    $deletion = "$ref[$j]";
	    $delstart = $j;
	    while (($j < $#ref) and ($alt[$j+1] eq "-")) {
		$j++;
		$deletion .= "$ref[$j]";
	    }
#	    print "deletion: $deletion\n";
	    $count[$delstart-$offset]{del}{$deletion} += $altcount;
	    $count[$delstart-$offset]{deltotal} += $altcount;
	} elsif ($alt[$j] eq ".") { # identity
	    # nothing to do
	} elsif ($alt[$j] ne $ref[$j] ) { # mismatch
	    $count[$j-$offset]{subtotal} += $altcount;
	} else { # shouldn't get here
	    print STDERR "Unhandled case: $ref[$j]\t$alt[$j]\n";
	}
	$j++;
    }
}

`echo -e "$reftrimid\t$totalcount" > $reftrimid.totalcount`;

# print out the results
#print "Reference Insertions Deletions Mismatches\n";
$indelfile = "$ifile.indel";
open INDELFILE, ">$indelfile" or die "cannot open $indelfile\n";
$datfile = "$ifile.dat";
open OFILE, ">$datfile" or die "cannot open $datfile\n";
print OFILE "Base\tType\tValue\n";
$labels = "";
$limits = "";
$vjust = "";
#for $j (0..$#refseq) {
$refseq[$#refseq+1] = "-"; # this allows for insertions at the end of the sequence
for $j (0..$#refseq) {
#    print "$j $refseq[$j] ";

#    for $i (0..2) {
#	if ($count[$i][$j] == 0) {
#	    $percent[$i] = 0;
#	} else {
#	    $percent[$i] = int($count[$i][$j] / $totalcount * 1000) / 10;
#	}
#    }	
    # calculate percent, set count to zero if undefined
    $count[$j]{intotal} = $count[$j]{intotal} // 0;
    $count[$j]{deltotal} = $count[$j]{deltotal} // 0;
    $count[$j]{subtotal} = $count[$j]{subtotal} // 0;
    $inpercent = &round1000($count[$j]{intotal} / $totalcount * 100);
    $delpercent = &round1000($count[$j]{deltotal} / $totalcount * 100);
    $subpercent = &round1000($count[$j]{subtotal} / $totalcount * 100);

    # special for horizon: calculate percent of each type of indel 
    foreach $insertion (keys %{ $count[$j]{in} }) {
	$pct = &round1000($count[$j]{in}{$insertion} / $totalcount * 100);
	print INDELFILE "$reftrimid\t$totalcount\t$j\tIN-$insertion\t$pct\n" if ($pct > 0.001);	
    }
    foreach $deletion (keys %{ $count[$j]{del} }) {
	$pct = &round1000($count[$j]{del}{$deletion} / $totalcount * 100);
	print INDELFILE "$reftrimid\t$totalcount\t$j\tDEL-$deletion\t$pct\n" if ($pct > 0.001);	
    }

    print OFILE "$j\tInsertions\t$inpercent\n";
    print OFILE "$j\tDeletions\t$delpercent\n";
    print OFILE "$j\tSubstitutions\t$subpercent\n";
    if (1) { #($j % 2 == 0) {
#	$xtics .= "\"$refseq[$j]\",";
	$labels .= "\"$j\" = \"$refseq[$j]\",";
	$limits .= "\"$j\",";
    } else {
	$labels .= "\"\\n$refseq[$j]\",";
    }
    $vjust .= ($j % 2) + 1 . ",";
}
$labels =~ s/,$//; # get rid of last comma
$limits =~ s/,$//; # get rid of last comma
$vjust =~ s/,$//; # get rid of last comma
close OFILE;

# run ggplot2 to generate plots
## Plotting ###
$ofile = $datfile;
$rfile = "$ifile.r";
$name = $ifile;
$name =~ s/\.demux\.fa//;

my $height = 480;
my $width = 2000;

### plot ###
#print "Generating plot\n";

open RFILE, ">$rfile" or die "Cannot open $rfile\n";
print RFILE qq(
  library(ggplot2);
  datat <- read.table("$ofile", header=T, colClasses=c("Base"="factor"));
  png(filename="$name.png", height = $height, width = $width);

  p <- ggplot(datat, aes(x=Base, y=Value))
  p + geom_bar(stat='identity', fill="red") + ylim(0,100) + xlab("Sequence") + ylab("Percent") + theme(legend.position="top") + theme(legend.title=element_blank()) + facet_wrap(~ Type, ncol = 1) + scale_x_discrete(labels=c($labels), limits=c($limits)) + theme(text = element_text(size = 20)) + theme(axis.text.x=element_text(size=7)) + labs(x="$name");
  dev.off();

  png(filename="$name-log.png", height = $height, width = $width);
  datat\$Value[datat\$Value < $sensitivitylimit] <- $sensitivitylimit
  p <- ggplot(datat, aes(x=Base)) 
  p + geom_linerange(aes(ymax=Value, ymin=$sensitivitylimit), size=2, color="red") + scale_y_log10(limits=c($sensitivitylimit,100), breaks=c(.1,.5,1,5,10,50,100)) + xlab("Sequence") + ylab("Percent") + theme(legend.position="top") + theme(legend.title=element_blank()) + facet_wrap(~ Type, ncol = 1) + scale_x_discrete(labels=c($labels), limits=c($limits)) + theme(axis.text.x=element_text(size=7)) + labs(x="$name");

  dev.off();
  #eof");

close RFILE;
$result = `R --no-restore --no-save --no-readline < $rfile &> $rfile.out`;
if ($?) {
    print STDERR "Error: R plotting issue: $result\n";
    $result = `cat $rfile.out`;
    print STDERR $result;
}

exit if ($verbose);

# clean up tmp files
`rm $altfa`;
`rm $reffa`;
`rm $rfile`;
`rm $datfile`;
`rm $rfile.out`;


# Round to neaerest hundreth
sub round1000 {
    return sprintf("%.3f", $_[0]);
#    return int($_[0] * 100 + 0.5) / 100;
}
