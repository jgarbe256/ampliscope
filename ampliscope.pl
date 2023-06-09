#!/usr/bin/perl -w

#####################################
# ampliscope.pl
# John Garbe
# June 2023
#####################################

=head1 DESCRIPTION

ampliscope.pl - Visualize variation in amplicon sequence data

=head1 SYNOPSIS

ampliscope.pl [--threads] [--margin 10] [--minreads 1] [--outputfolder folder] [--referencefasta amplicon.fasta] primers.txt R1.fastq R2.fastq

=head1 OPTIONS

This pipeline only handles paired-end reads. Fastq files may be gz compressed.

Options:

 --help : Print usage instructions and exit
 --verbose : Print more information while running
 --outputfolder string : output folder name (ampliscope-SAMPLENAME)
 --threads integer : Number of threads/cores to use (10)
 --margin integer : Keep sequences within +/- margin of the expected lenth (10)
 --minreads integer : Keep sequences with at least minreads reads (1)
 --revcomp : Reverse complement the second column of primer sequences
 --printlimit integer : # how many unique inserts to print out (10)
 --dimerlength integer : # sequences shorter than the primer lengths plus dimerlength is considered a primer dimer (10)
 --minidentity integer : Discard sequences with less than --minidentity identity to the reference sequence (default 0, set to 90 or 95 to remove off-target sequences)
 --referencefasta file : fasta file containing reference amplicon sequences. Sequence names must match amplicon names in primer file.

=cut

##################### Initialize ###############################

use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Pod::Usage;
use File::Temp qw( tempdir );

# set defaults
$args{threads} = 10;
$args{margin} = 10;
$args{minreads} = 1;
$args{verbose} = '';
$args{pearmin} = 10;
$args{printlimit} = 10;
$args{dimerlength} = 10;
$args{minidentity} = 0;
GetOptions("help" => \$args{help},
	   "verbose" => \$args{verbose},
	   "threads=i" => \$args{threads},
	   "minreads=i" => \$args{minreads},
	   "margin=i" => \$args{margin},
	   "outputfolder=s" => \$args{outputfolder},
	   "referencefasta=s" => \$args{referencefasta},
	   "revcomp" => \$args{revcomp},
	   "printlimit=i" => \$args{printlimit},
	   "dimerlength=i" => \$args{dimerlength},
	   "minidentity=i" => \$args{minidentity},
    ) or pod2usage;
pod2usage(-verbose => 99, -sections => [qw(DESCRIPTION|SYNOPSIS|OPTIONS)]) if ($args{help});

### Handle Parameters ###
if ($#ARGV != 2) {
    print "Expecting three commandline parameters\n";
    pod2usage;
}
$args{primerfile} = shift @ARGV;
$args{r1} = shift @ARGV;
$args{r2} = shift @ARGV if (@ARGV);
# get samplename from fastq file name
my ($samplename) = fileparse($args{r1});
$samplename =~ /(.*)_R1/;
$name = $1 // "";
if (! $name) {
    $samplename =~ /(.*)\.fast[qa]/;
    $name = $1 // "";
}
$samplename = $name || "sample";
print STDERR "processing sample $samplename\n";

die "Cannot find r1 file $args{r1}\n" unless (-e $args{r1});
$args{r1} = abs_path($args{r1});
if ($args{r2}) {
    die "Cannot find r2 file $args{r2}\n" unless (-e $args{r2});
    $args{r2} = abs_path($args{r2});
}
die "Cannot find primer file $args{primerfile}\n" unless (-e $args{primerfile});
$args{primerfile} = abs_path($args{primerfile});
$args{verbose} = "--verbose" if ($args{verbose});
die "Cannot find reference fasta file $args{referencefasta}\n" if ($args{referencefasta} and (! -e $args{referencefasta}));
$args{referencefasta} = abs_path($args{referencefasta}) if ($args{referencefasta});

$args{outputfolder} = "ampliscope-$samplename" unless ($args{outputfolder});
`mkdir $args{outputfolder}` unless (-e $args{outputfolder});
chdir $args{outputfolder} or die "Cannot chdir to $args{outputfolder}\n";
`mkdir scratch` unless (-e "scratch");
chdir "scratch" or die "Cannot chdir to scratch\n";

# check for dependencies
foreach $program ("mview", "parallel", "mafft", "pear", "needle", "R") {
    &depcheck($program);
}

### Merge reads with pear ###
print STDERR "Merging reads with pear:\n";
if ($args{r2}) {
    $length = "";
    if ($args{min} and $args{max}) {
	$length = "-n $args{min} -m $args{max}";
    } else {
	$length = "-n $args{pearmin}"; # this will help stitch primer dimers
    }
    $result = `pear -f $args{r1} -r $args{r2} -e -o $samplename $length -j $args{threads} > $samplename.pear.log`; # don't calculate empirical frequencies: calculation takes a long-ass time
    die "Pear failure: $result\n" if ($?);
    $fastqfile = "$samplename.assembled.fastq";
} else {
    print STDERR "single-end reads, skipping merging\n";
    $fastqfile = $args{r1};
}

# print out # and % of reads merged
$result = `grep "^Assembled reads \\." $samplename.pear.log`;
print $result;
$result =~ s/,//g;
$result =~ /: (\d+) \/ (\d+)/;
$mergedreads = $1;
$totalreads = $2;
#print "$mergedreads, $totalreads\n";

### Dedup with dedup.pl ###
print STDERR "Deduplicating reads with dedup.pl:\n";
$dedupedfasta = "$samplename.deduped.fasta";
`dedup.pl -m $args{minreads} $args{verbose} $fastqfile > $dedupedfasta`; # TODO: maybe do the minreads later on
die "dedup.pl failure\n" if ($?); 

### Demultiplex by primer ###
print STDERR "Demultiplexing sequences by primers:\n";
$revcomp = "";
$revcomp = "--revcomp" if ($args{revcomp});
$referencefasta = ($args{referencefasta}) ? "--referencefasta $args{referencefasta}" : "";
print "primerprocess.pl $args{verbose} $revcomp $referencefasta --dimerlength $args{dimerlength} --printlimit $args{printlimit} --margin $args{margin} $args{primerfile} $dedupedfasta > $samplename.txt\n" if ($args{verbose});
`primerprocess.pl $args{verbose} $revcomp $referencefasta --dimerlength $args{dimerlength} --printlimit $args{printlimit} --margin $args{margin} $args{primerfile} $dedupedfasta > $samplename.txt`;
if ($?) {
    die "primerprocess error\n";
}

### Generate gnuplots ###
print STDERR "Generating variation plots with align-parse.pl:\n";
$commandfile = "align-parse-commands.txt";
open OFILE, ">$commandfile" or die "cannot open $commandfile: $!\n";
@fas = `ls *.demux.fa`;
foreach $file (@fas) {
    chomp $file;
    $lines = `wc -l < $file`;
    chomp $lines;
    if ($lines == 0) {
	print STDERR "Skipping file with no sequences: $file\n";
	next;
    }
    print OFILE "align-parse.pl -i $args{minidentity} -v $file\n";
}
close OFILE;
`cat $commandfile | parallel -j $args{threads}`;
# TODO: add xargs alternative to parallel
# TODO: count how many reads were discarded

# align on all fasta files and generate alignment views
print STDERR "Generating html alignments:\n";
$commandfile = "html-alignments-commands.txt";
open OFILE, ">$commandfile" or die "cannot open $commandfile: $!\n";
@fas = `ls *.good.fa`;
foreach $file (@fas) {
    chomp $file;
    $lines = `wc -l < $file`;
    $name = $file;
    $name =~ s/\.good\.fa//;
    chomp $lines;
    if ($lines == 0) {
	print STDERR "Skipping file with no sequences: $file\n";
	next;
    }
    print OFILE "mafft --globalpair  $file > $file.aln 2> $file.aln.log && mview -in pearson -out new -html head -css on -coloring identity -ruler on -moltype dna $file.aln > $name.html\n";
}
close OFILE;
`which parallel`;
if ($?) {
    `cat $commandfile | xargs -L 1 -I CMD -P $args{threads} bash -c CMD`;
} else {
    `cat $commandfile | parallel -j $args{threads}`;
}

die "Error: mafft and/or mview failed on one or more amplicons. Take a look at *.aln.log files" if ($?);

### Generate html index ###
print STDERR "Generating html index...\n";
$ofile = "../$samplename.html";
open OFILE, ">$ofile" or die "cannot open output file $ofile: $!\n";
foreach $fa (sort hdsort @fas) {
    chomp $fa;
    $name = $fa;
    $name =~ s/\.good\.fa//;
    if (-e "$name.totalcount") {
	$totalcount = `cut -f 2 $name.totalcount | head`;
	chomp $totalcount;
	$png = "$name.png";
	print OFILE "<img src=\"png\\$png\"><br>\n";
	$html = "$name.html";
	print OFILE "<a href=\"html\\$html\">$name</a> $totalcount reads<br>\n";
    } else {
	print OFILE "$name 0 reads<br>\n";
    }
}

# log-scale version
$ofile = "../$samplename-log.html";
open OFILE, ">$ofile" or die "cannot open output file $ofile: $!\n";
foreach $fa (sort hdsort @fas) {
    chomp $fa;
    $name = $fa;
    $name =~ s/\.good\.fa$//;
    if (-e "$name.totalcount") {
	$totalcount = `cut -f 2 $name.totalcount | head`;
	chomp $totalcount;
	$png = "$name-log.png";
	print OFILE "<img src=\"png\\$png\"><br>\n";
	$html = "$name.html";
	print OFILE "<a href=\"html\\$html\">$name</a> $totalcount reads<br>\n";
    } else {
	print OFILE "$name 0 reads<br>\n";
    }
}

# print out a couple summary metrics
print "$totalreads total reads\n";
print "$mergedreads reads after merging with pear\n";
$mergedreadspct = &round10($mergedreads / $totalreads * 100);

# add badidentity column to demultiplex-stats.txt
@results = `cat *.badcount`;
foreach $result (@results) {
    chomp $result;
    ($amplicon, $value) = split /\t/, $result;
    $identityvalue{$amplicon} = $value;
}
open IFILE, "demultiplex-stats.txt" or die "cannot open stats file: $!\n";
open OFILE, ">ampliscope-stats.txt" or die "cannot open stats output file: $!\n";
$header = <IFILE>;
chomp $header;
print OFILE "$header\tBadIdentity\n";
while ($line = <IFILE>) {
    chomp $line;
    @line = split /\t/, $line;
    $amplicon = $line[0];
    $value = $identityvalue{$amplicon} // 0;
    print OFILE "$line\t$value\n";
}
close OFILE;
close IFILE;

# create plot of reads per amplicon
$name = "ampliconplot";
$ofile = "ampliscope-stats.txt";
$rfile = "$name.rmd";
print "Generating $name\n";

$title = "Amplicon Read Counts";
$text = "Number of reads per amplicon, after read stitching. Reads are considered dimers if there are fewer than $args{dimerlength}bp between the primer sequences. Reads are considered bad length if the stitched length is $args{margin}bp longer or shorter than the expected length. Reads are considered low identity if they have less than $args{minidentity} percent identity to the reference sequence.";
@metrics = ("DimerReads", "BadLengthReads", "BadIdentity", "PassedReads");
@metricnames = ("Dimers", "Bad Length", "Low Identity", "Passed Reads");
$traces = "";
for $metric (1..$#metrics) {
    $traces .= qq(%>% add_trace(y=~$metrics[$metric],  name = "$metricnames[$metric]", hovertext = paste0(datat\$sample,"<br>$metricnames[$metric] ",round(datat\$$metrics[$metric],2))) );
}


open RFILE, ">$rfile" or die "Cannot open $rfile\n";
print RFILE qq(
<b>$title</b>
<div id="$name">
$text
</div>

```{r $name}

library(plotly)

datat <- read.table("$ofile", header=T, colClasses=c("Amplicon"="factor"));

sampleorder <- list(title = "Sample", automargin = TRUE, categoryorder = "array", tickmode = "linear", tickfont = list(size = 10),
               categoryarray = sort(c(as.character(datat\$Amplicon))))

config(plot_ly(
  data=datat,
  x = ~Amplicon,
  y = ~$metrics[0],
 name = "$metricnames[0]",
 type = "bar",
hoverinfo="text",
hovertext = paste0(datat\$Amplicon,"<br>$metricnames[0] ",round(datat\$$metrics[0],2)) ) $traces
%>% layout(xaxis = sampleorder, dragmode='pan', legend = list(orientation = 'h'), barmode='stack', xaxis = list(title = "Sample"), yaxis = list(title = "Number of reads")), displaylogo = FALSE, modeBarButtonsToRemove = list('sendDataToCloud','toImage','autoScale2d','hoverClosestCartesian','hoverCompareCartesian','lasso2d','zoom2d','select2d','toggleSpikelines','pan2d'))
```
);
close RFILE;

# create pseudogel
@fas = `ls *.demux.fa`;
$ofile = "pseudogel.filelist";
open OFILE, ">$ofile" or die $!;
foreach $file (@fas) {
    chomp $file;
    $lines = `wc -l < $file`;
    $name = $file;
    $name =~ s/\.demux\.fa//;
    print OFILE "$file\t$name\n";
}
$result = `pseudogel.pl -l $ofile -d`;
print $result if ($args{verbose});

# Create report
my $date = `date`;
chomp $date;
$ampliconcount = `wc -l < demultiplex-stats.txt`;
chomp $ampliconcount;
$readtype = "Single-end";
if ($args{r2}) {
    $readtype = "Paired-end";
}

$reportfile = "report.rmd";
open $reportFH, ">$reportfile" or die $!;
print $reportFH qq(---
title: Ampliscope Analysis Report - $samplename
output: html_document
---

<style>
 .main-container {
   width: 95%;
   max-width: 800px;
  }
</style>

```{r setup, include=FALSE}
knitr::opts_chunk\$set(echo=FALSE, message=FALSE, error=TRUE, out.width="100%", out.height=300)
htmltools::tagList(rmarkdown::html_dependency_font_awesome())
```
<b>Sample name:</b> $samplename<br>
<b>Read type:</b> $readtype<br>
<b>Amplicons:</b> $ampliconcount<br>
<b>Total reads:</b> $totalreads<br>
<b>Reads after merging:</b> $mergedreads ($mergedreadspct%)<br>
<b>Report generated:</b> $date<br>

```{r test-main, child = 'pseudogelplot.rmd'}\n```

```{r test-main, child = 'ampliconplot.rmd'}\n```

[Variation Plots]($samplename.html)  
[Variation Plots log scale]($samplename.html)  
[Unmatched barcodes](unmatched-barcodes.txt)  

## Methods
Overlapping read-pairs were assembled with Pear, non-overlapping read-pairs and assembled reads shorter than $args{pearmin} were discarded. Reads were demultiplexed by amplicon based on exact matches to amplicon primer sequences. Reads with an assembled length $args{margin}bp longer or shorter than the expected length of the amplicon were discarded. For each amplicon the stitched reads were deduplicated into a list of unique sequences, with a count of the number of reads seen for each sequence. Sequences with a read count less than $args{minreads} were discarded. The most abundant sequence was selected to be the reference sequence for the amplicon.  Needle (Emboss package) was used to generate optimal global sequence alignments between each amplicon sequence and the amplicon reference sequence. The numbers of insertions, deletions, and substitutions at each base of the reference amplicon sequence were counted. Alignments of the $args{printlimit} most common amplicon reads were visualized using MView.

);

# render report into html
$result = `Rscript -e \"library('rmarkdown'); render('$reportfile', output_format='html_document', output_file='../report.html')\" 2>&1`; 
if ($?) {
    print STDERR "Error generating report: $result\n";
} else {
    print STDERR $result if ($args{verbose});
    print STDERR "Report file report.html created\n";
}

### Clean up ###
# move keepers out of scratch folder
`mkdir ../png` unless (-e "../png");
`mv *.png ../png`;

`mkdir ../html` unless (-e "../html");;
`mv *.html ../html`;

`mkdir ../fa` unless (-e "../fa");
`mv *.fa ../fa`;

#`mkdir ../align` unless (-e "../align"); # these can stay in scratch
#`mv *.aln ../align`;

`mv unmatched-barcodes.txt ../`;

print STDERR "Finished. Have a nice day!\n";
exit;

################################ Helper subs ###############################

### sort function for horizon discovery
sub hdsort {
    my ($mya) = split /\./, $a;
    my ($myb) = split /\./, $b;
    return $mya cmp $myb;
}

# Round to neaerest tenth
sub round10 {
    return sprintf("%.1f", $_[0]);
}

# Round to neaerest hundreth
sub round100 {
    return sprintf("%.2f", $_[0]);
}

sub mean {
    my($data) = @_;
    if (not @$data) {
	die("Empty array\n");
    }
    my $total = 0;
    foreach (@$data) {
	$total += $_;
    }
    my $average = $total / @$data;
    return $average;
}

sub stdev {
    my($data) = @_;
    if(@$data == 1) {
	return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach(@$data) {
	$sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / (@$data-1)) ** 0.5;
    return $std;
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}

sub depcheck {
    my $program = shift @_;
    my $result = `which $program 2>&1`;
    die "Error, required program $program not found, please refer to ampliscope installation instructions.\n $result\n" if ($?);
}
