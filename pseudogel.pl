#!/usr/bin/perl -w

#######################################################################
# pseudogel.pl
# John Garbe
# April 2023
#
#######################################################################

=head1 NAME

pseudogel.pl - Generate a pseudogel of library fragment lengths

=head1 SYNOPSIS

pseudogel.pl -l filelist

=head1 DESCRIPTION

Generate a plot showing distribution of fragment lengths per sample

 -l file : a file with a list of fastq or fastq files, one per line. A second tab-delimited column contains the sample names
 -d : The reads have been deduplicated, pull read counts from the read name
 -h : Print usage instructions and exit
 -v : Print more information while running (verbose)

=head1 EXAMPLE

Run the script::

    $ pseudogel.pl -l filelist.txt

=cut

############################# Main  ###############################

use Cwd 'abs_path';
use Getopt::Std;
use Term::ANSIColor;
use Pod::Usage;
use File::Temp qw/ tempdir /;

our ($opt_h, $opt_v, $opt_l, $opt_f, $opt_d);

$usage = "USAGE: pseudogel.pl -l filelist\n";
die $usage unless getopts('hvl:fd');
die $usage unless ($opt_l);
pod2usage(q(-verbose) => 3) if ($opt_h);

# Range of the gel, and granularity - maybe don't need this
#$lengthmin = 0;
#$lengthmax = 1000;
#$bin = 10;

### parameters
$filelist = abs_path($opt_l) unless (!$opt_l);
die "A filelist must be provided using the -l option\n" unless ($filelist);
die "Filelist $filelist not found\n" if ($filelist && (!-e $filelist));
my $verbose = $opt_v // 0; # print out more stuff

$numsamples = 0;

### Run through file list ###
open $fl, "$filelist" or die "Cannot open filelist $filelist: $!\n";
while ($fline = <$fl>) {
    chomp $fline;
    ($file, $sample) = split /\t/, $fline;
    if (! -e $file) {
	print "File $file not found, skipping\n";
	next;
    }
    print "Processing file $file\n" if ($verbose);
    $numsamples++;

    # open the file, compressed or not
    if ($file =~ /\.gz$/) {
	open IFILE, "gunzip -c $file |" or die "Can’t open gzipped input file $file: $!\n";
    } else {
	open IFILE, "$file" or die "Cannot open input file $file: $!\n";
    }

    # get size info from input file
    $totalsequences = 0;
    open IFILE, $file or die "Cannot open $file: $!\n";
    if (($file =~ /\.fastq/) or ($file =~ /\.fq/)) { # fastq format
	while ($line = <IFILE>) {
	    $seq = <IFILE>; # seq
	    $line = <IFILE>; # plus
	    $line = <IFILE>; # qual
	    chomp $seq;
	    $length = length($seq);
	    $data{$sample}{$length}++;
	    $totalsequences++;
	}
	
    } else { # assume fasta format
	$seq = "";
	$readcount = 1;
	while ($line = <IFILE>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		if ($seq) {
		    $length = length($seq);
		    $data{$sample}{$length} += $readcount;
		    $totalsequences += $readcount;
		    $seq = "";
		} 
		if ($opt_d) { # parse copy number from sequence name
		    ($readcount) = $line =~ /:(\d+)$/;
		}
	    } else {
		$seq .= $line;
	    }
	}
	# handle last sequence
	$length = length($seq);
	$data{$sample}{$length} += $readcount;
	$totalsequences += $readcount;
    }
    close IFILE;

    # calculate stats
    $stats{$sample}{"totalsequences"} = $totalsequences;
}
close $fl;

# print out stats
print "STATS Start\n";
foreach $sample (keys %stats) {
    foreach $key (keys %{$stats{$sample}}) {
	print "$sample\t$key\t$stats{$sample}{$key}\n";
    }
}
print "STATS End\n";

# DEBUG: print out data hash, and convert from read count to percent
foreach $sample (keys %data) {
    print STDERR "\n$sample\n";
    foreach $length (sort {$a <=> $b} keys %{$data{$sample}}) {
#	print "$length\t$data{$sample}{$length}\n";
	if ($stats{$sample}{totalsequences}) {
	    $percent = int($data{$sample}{$length} / $stats{$sample}{totalsequences} * 100);
	} else {
	    $percent = 0;
	}
	$percent{$sample}{$length} = $percent;
	print STDERR "$length\t$percent\n";
    }
}

# generate cycleplot plots - not ready to test this yet
$bins = int(100/($numsamples/10));
$title = "Library Insert Size Pseudogel";
&cycleplot("pseudogelplot", \%percent, $title, $bins);

exit;

####################### helper subs #######################

sub cycleplot {
    my ($name, $data, $title, $bins) = @_;
    
    my (%smalldata, $cycles, @samples, @cycles, %binsize, %color);
    
    ### parameters
    
    $ofile = "$name.dat";
        
    @samples = sort keys %{$data};
#    @cycles = sort {$a<=>$b} keys %{$data->{$samples[0]}};
#    $cycles = $cycles[$#cycles];
    $cycles = 100; # hardcode the cycles to 1kb for now
    $maxbins = int($cycles / 10); # put a limit on the number of bins
    $bins = $maxbins if ($bins > $maxbins);
    $bins = 5 if ($bins < 5);
    $bins = 100; # DEBUG
    print STDERR "$bins bins\n";
    
    # cluster into bins with ckmeans
    # print data out to tmp file, run ckmeans on it
    $ofile = "pseudogel.dat";
    open $OFILE, ">$ofile" or die "cannot open $ofile: $!\n";
    foreach $sample (@samples) {
	foreach $cycle (1..$cycles) {
	    $value = $data->{$sample}{$cycle} // 0;
	    print $OFILE "$value ";
	    $total{$sample} += $value;
	}
	print $OFILE "\n";
    }
    close OFILE;

    print "ckmeans.r $ofile $bins\n";
    @results = `ckmeans.r $ofile $bins`;
#    print @results;

    chomp @results;
    # process results
    foreach $sample (@samples) {
	$clusters = shift @results;
	$centers = shift @results;
	@clusters = split /\t/, $clusters;
	@centers = split /\t/, $centers;
	# calculate size of each bin
	foreach $cluster (@clusters) {
	    $binsize{$sample}{$cluster}++;
	}
	# store value (center) for each bin
	for $i (0..$bins-1) {
	    $bin = $i+1;
	    $smalldata{$sample}{$bin} = $centers[$i] // 0;
	    $binsize{$sample}{$bin} = 0 unless ($binsize{$sample}{$bin});
	}
    }
 
    # print out data
    $ofile = "$name.dat";
    open OFILE, ">$ofile" or die "cannot open temporary file $ofile for writing: $!\n"; 
    print OFILE "sample";
    for $bin (1..$bins) {
	print OFILE "\tcluster$bin";
    }
    print OFILE "\ttotal\n";
    foreach $sample (@samples) {
	print OFILE "$sample";
	for $bin (1..$bins) {
	    print OFILE "\t$binsize{$sample}{$bin}";
	}
	print OFILE "\t$total{$sample}\n";
    }
    close OFILE;


    $colors = qq(
colfunc<-colorRampPalette(c('black', 'white'))
mycolors <- rev(colfunc(n=100+1))
);
#    $altcolors = qq(
#colfunc<-colorRampPalette(c('darkred','indianred', 'yellow', 'gold', 'darkgreen', 'green'))
#mycolors <- colfunc(n=41+1)
#); # another option, may be useful
    $sort = "data2 <- datat[order(datat\$total),]";
#   $sort = "data2 <- datat[rev(order(datat\$total)),]";	
    $preunit = "";
    $postunit = "% Adapter";
    $text = "This plot shows the library insert fragment lengths for each sample.";
    $ordername = "Insert length";

    # generate colors and text
    for $bin (1..$bins) {
	$color{$bin} = "marker = list(color = c(";
	$text{$bin} = "c(";
	foreach $sample (@samples) {
	    print "undefined: $sample $bin\n" if (!defined($smalldata{$sample}{$bin}));
	    $color{$bin} .= "mycolors[$smalldata{$sample}{$bin}+1],"; # colors are 1-based
	    $value = round10($smalldata{$sample}{$bin});
#	    $maxadapter = $maxadapter{$sample};
#	    $maxadapter =~ s/'/\\'/g;
#	    $aname = ($name =~ /fastqadapterplot/) ? "<br>$maxadapter" : "";
	    $aname = "aname";
	    $text{$bin} .= "'$sample<br>$preunit$value$postunit$aname',";

	}
	chop $color{$bin}; # take off last ,
	$color{$bin} .= "))";
	chop $text{$bin}; # take off last ,
	$text{$bin} .= ")";
    }	

    $updatemenues = "updatemenus = updatemenus,";
    $updatemenues = ""; #disable sort buttons

    # plot the data
    $traces = "";
    for $bin (2..$bins) {
	$traces .= " %>% add_trace(y=~cluster$bin, $color{$bin}, hovertext = $text{$bin}) ";
    }

    $rfile = "$name.rmd";
    open RFILE, ">$rfile" or die "Cannot open $rfile\n";
    print RFILE qq(
<div class="plottitle">$title
<a href="#$name" class="icon" data-toggle="collapse"><i class="fa fa-question-circle" aria-hidden="true" style="color:grey;font-size:75%;"></i></a></div>
<div id="$name" class="collapse">
$text
</div>

```{r $name}

library(plotly)

datat <- read.table("$ofile", header=T, colClasses=c("sample"="factor"));

$colors

sampleorder <- list(title = "Sample", automargin = TRUE, categoryorder = "array", tickmode = "linear", tickfont = list(size = 10), 
               categoryarray = sort(c(as.character(datat\$sample))))

$sort
valueorder <- list(title="Sample", automargin = TRUE, categoryorder = "array", tickmode = "linear", tickfont = list(size = 10),
              categoryarray = c(as.character(data2\$sample)))

updatemenus <- list(
  list(
    active = 0,
    font = list(size = 10),
    type = 'buttons',
    xanchor = 'right',
    buttons = list(
      
      list(
        label = "Sort by<br>Sample",
        method = "relayout",
        args = list(list(xaxis = sampleorder))),
      
      list(
        label = "Sort by<br>Value",
        method = "relayout",
        args = list(list(xaxis = valueorder)))
    )
  )
)

config(plot_ly(
  data=datat,
  x = ~sample,
  y = ~cluster1,
 type = 'bar', showlegend=FALSE, $color{1},
hoverinfo="text",
hovertext = $text{1} )
$traces
%>% layout(xaxis = sampleorder, $updatemenues xaxis = list(title='Sample'), dragmode='pan', showlegend = FALSE, barmode='stack', yaxis = list(title = "Insert Size (bp)", fixedrange = TRUE, showgrid = FALSE)), displaylogo = FALSE, modeBarButtonsToRemove = list('sendDataToCloud','toImage','autoScale2d','hoverClosestCartesian','hoverCompareCartesian','lasso2d','zoom2d','select2d','toggleSpikelines','pan2d','collaborate'))
```
);

    # Print out samples ordered by duplication - todo: calculate median length
    @order = sort {$total{$b} <=> $total{$a}} keys %total;
    $order = "";
    foreach $sample (@order) {
	$order .= "\"$sample\",";
    }
    chop $order; # remove last ,
    print "ORDER\t$ordername\t$order\n";

}

########## subs ###############

# Round to neaerest tenth
sub round10 {
    return sprintf("%.1f", $_[0]);
}

# Round to neaerest hundreth
sub round100 {
    return sprintf("%.2f", $_[0]);
}

sub round_to_bin {
    my ($num, $bin) = @_;
    my $rounded_num = int($num/$bin + 0.5) * $bin;
    return $rounded_num;
}
