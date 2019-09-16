Small RNA Seq - miR-210 Practical
================
Jack Monahan, Anton Enright & Katy Brown

Feb 2018

-   [Analysis of an sRNA-Seq data for miR-210 over-expression](#analysis-of-an-srna-seq-data-for-mir-210-over-expression)
    -   [Experiment Overview](#experiment-overview)
    -   [Raw Data](#raw-data)
    -   [Count loading & Normalisation](#count-loading-normalisation)
    -   [Post Normalisation QC](#post-normalisation-qc)
    -   [Variance Stabilising the Counts](#variance-stabilising-the-counts)
    -   [Analysis of treatment effects](#analysis-of-treatment-effects)
    -   [Statistical Analysis - Differential Expression](#statistical-analysis---differential-expression)
    -   [Final Results](#final-results)

Analysis of an sRNA-Seq data for miR-210 over-expression
========================================================

Experiment Overview
-------------------

This set of sequencing data is from libraries prepared in a previous course. We have matched microRNA and mRNA samples from   MCF7 Human Breast Cancer Cells. The experiment here was to overexpress miR-210 to explore the effect on mRNAs. Control samples used a non-mouse microRNA from c.elegans (cel-miR-67) TCACAACCTCCTAGAAA

| filename                       | group | samplename | treatment | fullname          | 3p-ad               |
|--------------------------------|-------|------------|-----------|-------------------|---------------------|
| 20127\_1\#1.converted.fastq.gz | 3     | 3A         | miR210    | miRNA03\_MIR\_GR3 | TGGAATTCTCGGGTGCCAA |
| 20127\_1\#2.converted.fastq.gz | 4     | 4A         | miR210    | miRNA04\_MIR\_GR4 | TGGAATTCTCGGGTGCCAA |
| 20127\_1\#3.converted.fastq.gz | 6     | 3B         | miR210    | miRNA14\_MIR\_GR6 | TGGAATTCTCGGGTGCCAA |
| 20127\_1\#4.converted.fastq.gz | 3     | 7A         | Scr       | miRNA07\_SCR\_GR3 | TGGAATTCTCGGGTGCCAA |
| 20127\_1\#5.converted.fastq.gz | 4     | 8A         | Scr       | miRNA08\_SCR\_GR4 | TGGAATTCTCGGGTGCCAA |
| 20127\_1\#6.converted.fastq.gz | 6     | 6B         | Scr       | miRNA16\_SCR\_GR6 | TGGAATTCTCGGGTGCCAA |

Raw Data
--------

The data should be preinstalled. Otherwise download the [raw data](http://wwwdev.ebi.ac.uk/enright-srv/courses/rna_cambridge_2017/course_counts/data) and put the files in a directory called 'mir210_smallRNA' in your main course_data directory.

<br>
Launch RStudio

``` r
setwd('~/Desktop/Course_Materials/small_RNASeq')
library(RColorBrewer)
library(gplots)
library(DESeq2)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
```

We can now load the count data

``` r
mircounts <- read.table("mircounts.txt",header=TRUE,row.names=1)
```

As well as the pdata, which contains information on each sample.

``` r
pdata <- read.table("pdata.txt",header=TRUE,row.names=1,sep="\t")
colnames(mircounts)=rownames(pdata)

groups=as.factor(pdata$group)
conds=as.factor(pdata$treatment)
```

Count loading & Normalisation
-----------------------------

We are now ready to create a DESeq object from the counts table. We are only looking at one condition, tretment of the samples with estradiol. We use the design formula “~ treatment”. We are going to estimate coefficients for the treatment conditions (Scrambled Control versus miR210 Overexpression)

``` r
coldata = as.data.frame(t(t(conds)))
rownames(coldata)=colnames(mircounts)
colnames(coldata)='treatment'

dds <- DESeqDataSetFromMatrix(countData = mircounts, colData = coldata, design = ~ treatment)
```

We are ready to normalise the data, but first we should look at the number of sequenced reads per sample.

``` r
cond_colours = brewer.pal(length(unique(conds)),"Accent")[as.factor(conds)]
names(cond_colours)=conds

group_colours = brewer.pal(3,"Accent")[as.factor(pdata$group)]
names(group_colours)=pdata$group
```

