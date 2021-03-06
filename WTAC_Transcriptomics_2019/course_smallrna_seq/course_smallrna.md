Small RNA Seq - Course Data
================
Anton Enright & Jack Monahan
'26 June, 2018'

-   [Analysis of smallRNA Seq datasets](#analysis-of-smallrna-seq-datasets)
    -   [Experimental Design](#experimental-design)
    -   [Setup](#setup)
    -   [SmallRNA Read Quality Control](#smallrna-read-quality-control)
    -   [Mapping Cleaned Reads to MicroRNAs](#mapping-cleaned-reads-to-micrornas)
-   [Analysis of smallRNA count data](#analysis-of-smallrna-count-data)
    -   [Experiment Setup](#experiment-setup)
    -   [Preparation](#preparation)
    -   [Main Analysis](#main-analysis)
    -   [Count Loading](#count-loading)
    -   [Count Preparation & Normalisation](#count-preparation-normalisation)

Analysis of smallRNA Seq datasets
=================================

In this practical we will be taking the 12 samples you generated libraries for on the course.

We will trim off the 3' Sequencing adapters and then you will construct a protocol for their analysis

The fastq files, pdata file and mircounts file are now installed on your desktops.

Experimental Design
-------------------

This is the sample description file used for the analyses below.

| **samplename** | **filename**               | **3p-ad**        | **genotype** | **group** |
|----------------|----------------------------|------------------|--------------|-----------|
| samplename     | filename                   | 3p-ad            | genotype     | group     |
| wt1            | 26151\_merged1.cram.fq.gz  | AGATCGGAAGAGCACA | wt           | g1        |
| wt2            | 26151\_merged2.cram.fq.gz  | AGATCGGAAGAGCACA | wt           | g4        |
| wt3            | 26151\_merged3.cram.fq.gz  | AGATCGGAAGAGCACA | wt           | g6        |
| nrep\_scram1   | 26151\_merged4.cram.fq.gz  | AGATCGGAAGAGCACA | nrep\_scram  | g1        |
| nrep\_scram2   | 26151\_merged5.cram.fq.gz  | AGATCGGAAGAGCACA | nrep\_scram  | g2        |
| nrep\_scram3   | 26151\_merged6.cram.fq.gz  | AGATCGGAAGAGCACA | nrep\_scram  | g3        |
| nrepm29a\_ho1  | 26151\_merged7.cram.fq.gz  | AGATCGGAAGAGCACA | nrepm29a\_ho | g2        |
| nrepm29a\_ho2  | 26151\_merged8.cram.fq.gz  | AGATCGGAAGAGCACA | nrepm29a\_ho | g4        |
| nrepm29a\_ho3  | 26151\_merged9.cram.fq.gz  | AGATCGGAAGAGCACA | nrepm29a\_ho | g5        |
| cyrano\_ho1    | 26151\_merged10.cram.fq.gz | AGATCGGAAGAGCACA | cyrano\_ho   | g3        |
| cyrano\_ho2    | 26151\_merged11.cram.fq.gz | AGATCGGAAGAGCACA | cyrano\_ho   | g5        |
| cyrano\_ho3    | 26151\_merged12.cram.fq.gz | AGATCGGAAGAGCACA | cyrano\_ho   | g6        |

Setup
-----

Lets load some libraries to get started

``` r
library(Reaper)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(ggplot2)
```

Now we will set our working directory to where the solexa FASTQ files (zipped) are stored

``` r
setwd("~/Desktop/course_data/wtac_smallrna_2018")
list.files()
```

               ## [1] "go.tar.gz"     "mircounts.txt" "pdata.txt"     "reaper.pdf"

Hopefully, you will see a compressed FASTQ txt file for each of the 4 lanes

It is important that we also load information for reaper that tells it the following:

-   Which FASTQ files are present.
-   Which Barcode sequences correspond to which sample names.
-   What 5' and 3' sequencing adapters were used in library generation.

You should see a sample table loaded into R:

``` r
pdata <- read.table("pdata.txt",header=TRUE,check.names=FALSE)
pdata
```

               ##      samplename                  filename            3p-ad    genotype
               ## 1           wt1  26151_merged1.cram.fq.gz AGATCGGAAGAGCACA          wt
               ## 2           wt2  26151_merged2.cram.fq.gz AGATCGGAAGAGCACA          wt
               ## 3           wt3  26151_merged3.cram.fq.gz AGATCGGAAGAGCACA          wt
               ## 4   nrep_scram1  26151_merged4.cram.fq.gz AGATCGGAAGAGCACA  nrep_scram
               ## 5   nrep_scram2  26151_merged5.cram.fq.gz AGATCGGAAGAGCACA  nrep_scram
               ## 6   nrep_scram3  26151_merged6.cram.fq.gz AGATCGGAAGAGCACA  nrep_scram
               ## 7  nrepm29a_ho1  26151_merged7.cram.fq.gz AGATCGGAAGAGCACA nrepm29a_ho
               ## 8  nrepm29a_ho2  26151_merged8.cram.fq.gz AGATCGGAAGAGCACA nrepm29a_ho
               ## 9  nrepm29a_ho3  26151_merged9.cram.fq.gz AGATCGGAAGAGCACA nrepm29a_ho
               ## 10   cyrano_ho1 26151_merged10.cram.fq.gz AGATCGGAAGAGCACA   cyrano_ho
               ## 11   cyrano_ho2 26151_merged11.cram.fq.gz AGATCGGAAGAGCACA   cyrano_ho
               ## 12   cyrano_ho3 26151_merged12.cram.fq.gz AGATCGGAAGAGCACA   cyrano_ho
               ##    group
               ## 1     g1
               ## 2     g4
               ## 3     g6
               ## 4     g1
               ## 5     g2
               ## 6     g3
               ## 7     g2
               ## 8     g4
               ## 9     g5
               ## 10    g3
               ## 11    g5
               ## 12    g6

Next we will start the Reaper algorithm. It will perform the following functions on all the lanes we have provided:

-   Splitting up reads according to provided barcodes
-   Detection of 3' or 5' adapter contamination using Smith-Waterman local alignment
-   Detection of low-complexity sequences, such as Poly-As or Poly-Ns
-   Quality score thresholding and trimming if required
-   Collapsing of reads according to depth in a summary FASTA result file per barcode per lane.
-   Generation of Quality Control plots for assessing sequencing quality

Reaper is started by passing it our samples table and telling it which "mode" to run in, in this case the mode is set to: **no-bc**.

    reaper(pdata,"no-bc");

Cleaning many millions of reads will take some time, the method processes around 2M-4M reads per minute and will look something like this:

    [1] "Starting Reaper for file: A638S1.R1.fastq.gz"
                   fastq                 geom                 meta             basename 
    "A638S1.R1.fastq.gz"              "no-bc"      "metadata1.txt"                  "1" 
                      do 
                 "10000" 
    Passing to reaper: dummy-internal --R -fastq A638S1.R1.fastq.gz -geom no-bc -meta metadata1.txt -basename 1 -do 10000

    ---
    mRpm   million reads per minute
    mNpm   million nucleotides per minute
    mCps   million alignment cells per second
    lint   total removed reads (per 10K), sum of columns to the left
    25K reads per dot, 1M reads per line  seconds  mr mRpm mNpm mCps {error qc  low  len  NNN tabu nobc cflr  cfl lint   OK} per 10K

Reaper is designed to be fast and memory-efficient so it should run on any machine with 500MB of RAM or more. The time taken to complete the run depends on how fast the processors in your machine are.

Let's take a look at the Quality Control Metrics generated

    reaperQC(pdata)

SmallRNA Read Quality Control
-----------------------------

Lets also make a nice PDF of the results and explore the QC metrics for the data

    pdf("reaper.pdf",width=12)
    reaperQC(pdata)
    dev.off()

You should get one plot back for each lane processed.

Mapping Cleaned Reads to MicroRNAs
----------------------------------

We can now use ChimiRa to map cleaned and filtered reads against miRBase known mature miRNA sequences. This server will take each cleaned unique read and compare it to precursor sequences downloaded from miRBase. In the case of mismatches, a read will still be assigned if it has less than two mismatches assuming both sequences are each others best hit. The system supports compressed (GZ or ZIP) files as well as FASTQ text files.

Click [here](http://www.ebi.ac.uk/research/enright/software/chimira) To use Chimira

For this web-server, choose each of the **.clean.uniq** fastq files that were produced by reaper. Make sure you choose **Human** as the reference species.

We now have raw counts of reads on microRNAs ready for QC and differential analysis.

Analysis of smallRNA count data
===============================

Experiment Setup
----------------

All data were pre-processed using *minion* to identify and check adapters, *reaper* to trim adapter sequences followed by *tally* to deduplicate reads while maintaining depth information. Subsequent to this all reads passed through the *ChimiRa* pipeline against all miRBase (Release 22) precursor sequences for Mouse. Reads were summed across paired end sequences for the same read pair. Finally reads are loaded into R for final analysis.

For this analysis, we will get you started.....

Preparation
-----------

Lets prepare some colour schemes for heatmaps and sample colours

``` r
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")(100)
```

Main Analysis
-------------

Count Loading
-------------

We can now load all the count data

``` r
mircounts <- read.table("mircounts.txt",header=TRUE,row.names=1)
mircounts=mircounts[-nrow(mircounts),]
```

Lets set up the count matrix column names and make a list of genotype conditions

``` r
colnames(mircounts)=rownames(pdata)
conds=as.factor(as.character(pdata$genotype))
conds
```

               ##  [1] wt          wt          wt          nrep_scram  nrep_scram 
               ##  [6] nrep_scram  nrepm29a_ho nrepm29a_ho nrepm29a_ho cyrano_ho  
               ## [11] cyrano_ho   cyrano_ho  
               ## Levels: cyrano_ho nrep_scram nrepm29a_ho wt

Count Preparation & Normalisation
---------------------------------

We are now ready to create a DESeq object from the counts table.

``` r
#Lets Load the Counts First
coldata = as.data.frame(t(t(conds)))
rownames(coldata)=colnames(mircounts)
colnames(coldata)='treatment'
dds <- DESeqDataSetFromMatrix(countData = mircounts, colData = coldata, design = ~ treatment)
```

We are ready to normalise the data, but first we should look at the number of sequenced reads per sample. There are some stark differences across the samples.

    cond_colours = brewer.pal(length(levels(conds)),"Set2")[conds]
    names(cond_colours)=pdata$genotype

    barplot(apply(mircounts,2,sum), las=2,col=cond_colours,main="Pre Normalised Counts",cex.names=0.4)
    legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
