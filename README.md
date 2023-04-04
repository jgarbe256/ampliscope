# ampliscope: Visualization of variation in amplicon sequencing data

------

The ampliscope packege visualizes variation in amplicon sequencing data. It takes raw sequencing reads in fastq format and a list of amplicon primers and generates plots and alignments showing the sequence variation present in the data for each amplicon.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Running Ampliscope](#running_ampliscope)
- [Ampliscope Output](#output)
- [Frequently Asked Questions](#FAQ)

## Overview<a name="overview"></a>

Ampliscope requires paired-end sequence data whose reads overlap at least 10bp to enable merging of the read-pairs. Input fastq files can be uncompressed or gz compressed.

Ampliscope implements the following analysis steps:

1: Overlapping read-pairs were assembled with Pear, non-overlapping read-pairs and assembled reads shorter than --minlength are discarded. 
2: Reads are demultiplexed by amplicon based on exact matches to amplicon primer sequences. Reads with an assembled length 10bp longer or shorter than the expected length of the amplicon were discarded. 
3: For each amplicon the stitched reads are deduplicated into a list of unique sequences, with a count of the number of reads seen for each sequence. Sequences with a read count less than 1 were discarded. The most abundant sequence is selected to be the reference sequence for the amplicon. 
4: The 10 most common amplicon reads are aligned to the reference amplicon sequence with the mafft multiple sequence aligner, and a visualization of the alignments is generated with MView
5: Needle (Emboss package) is used to generate optimal global sequence alignments between each amplicon sequence and the amplicon reference sequence. 
6: The needle alignments are parsed and the numbers of insertions, deletions, and substitutions at each base of the reference amplicon sequence are counted and plots summarizing the counts are generated.

## Installation<a name="installation"></a>

Install the following dependencies:
* Perl
* [`emboss`](<https://emboss.sourceforge.net/download/>) Alignment tool
* [`R`](<https://r-project.org>) The libraries plotly and rmarkdown are required
* [`mafft`](<https://mafft.cbrc.jp/alignment/software/>) Multiple sequence alignment tool
* [`pear`](<https://github.com/tseemann/PEAR>) Paired-end read meger
* [`mview`](<https://desmid.github.io/mview/>) Multiple sequence alignment visualization
* [`pandoc`](<https://pandoc.org>) Document converter
* [`GNU Parallel`](<https://www.gnu.org/software/parallel/>) Parallel job execution

Obtain a copy of the ampliscope package source code. You can either download and unzip the latest source from the github [release page](https://github.com/aryeelab/ampliscope/releases), or you use git to clone the repository by running `git clone --recursive https://github.com/jgarbe/ampliscope.git`


### Running Ampliscope<a name="running_ampliscope"></a>

Create a primer file listing the name of each amplicon, the forward and reverse primer sequences, and (optionally) the expected length of the amplicon. The file should be tab-delimited plain text:

>amplicon_01	ACAACGTTAGCCTGTT GTTGATATCCCACCCGAA	47
>amplicon_17	CCAAAAACAACAGTCA ATGGTGCCATTCTCCTT	115
>ldlr_1a	TTATCTGCTTGCTTCTGC	 ACTCCTGCAGGTCACTG	163
>ldlr1b	TTGTTAGGATGGTGGA	CAGGGCCTTTCCTCGC	81
>chr17-156890-v1	AGCCGGGACCACCT	TGGAGGTGAGGGAGAGG	240
>chr17-156890-alt	GGCCCGACTTGCAACTA	CTACCGGAGACGTGTCA	173

Run ampliscope:

```
ampliscope.pl --threads 10 primers.txt sample_R1.fastq.gz sample_R2.fastq.gz
```

Ampliscope options:

    This pipeline only handles paired-end reads. Fastq files may be gz
    compressed.

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

### Output<a name="output"></a>

The ampliscope output folder will contain these items:

- `report.html`: The main report, with link to all html pages and png images generated by ampliscope.
- `scratch`: Temporary files generated during the analysis.
- `png`: Plots for each amplicon summarizing insertions, deletions, and substitutions at each position of the amplicon.
- `html`: A web page for each amplicon visualizing the alignment of the most common sequences agains the reference.
- `fa`: A fasta file for each amplicon listing the most common sequences.
- `unmatched-barcodes.txt`: A list of the most common first 20bp of the R1 and R2 reads that didn't match any primer sequences listed in the primer file.

## Frequently Asked Questions<a name="FAQ"></a>

### Can I analyze data with UMIs

No, the pipeline does not support UMIs.

### What if my data is already demultiplexed and I have many fastq files, each containing reads for one amplicon?

Ampliscope isn't designed to directly handle this situation. You can combine all of your samples into a single pair of fastq files and then run Ampliscope.