Mini-mapping Practical
================
Anton Enright, Jack Monahan, Katy Brown

Basic Mapping Practical
-------

For any RNASeq experiment there must be a _mapping_ and a _quatification_ to turn mapped reads into counts on genomic features.

## Pre-requisites

For most mapping tools you will need:
* A genome (usually in FASTA format) - [Ensembl](http://www.ensembl.org/info/data/ftp/index.html) is a great place for this.
* An annotation file (usually in GTF format) - again, [Ensembl](http://www.ensembl.org/info/data/ftp/index.html) is ideal for this.
  * If you are using mouse or human data you could also look at [Gencode.](https://www.gencodegenes.org)
* An indexed version of the genome - the command needed varies by method, but the index needs to be built only once.
* A mapping tool (in this case HISAT2)
* A quantification tool (in this case htseq-count) to turn mapped reads to counts on features.


## Raw Data
The data should be preinstalled. Otherwise download the [raw data](http://wwwdev.ebi.ac.uk/enright-srv/courses/rna_cambridge_2017/mapping/data) and put the files in a directory called 'mir210_mapping' in your main course_data directory.

## Initial Setup

We have made a simple test with only Chromosome 22 and 400k reads on your machines

We will open a terminal and change to the practical folder on your desktop. 
We will also add the paths to the programs we will be using **hisat2** and **htseq-count**.

```
cd ~/Desktop/course_data/mapping_practical/

```

## Commonly used assembly and transcriptome analysis tools

| Method | Known Transcriptome | de novo Assembly | Ab initio Assembly |
|--------|---------------------|------------------|--------------------|
| Cufflinks | Yes | Yes | No |
| Trinity | No | No | Yes |
| Scripture | Yes | Yes | No |
| TransAbyss | No | No | Yes |
| Velvet | No | Yes |
| Rum | Yes | Yes | No |

Many of these tools require you to perform your own mapping using a mapping algorithm

## Commonly used mapping tools

*   Bowtie2
*   GMAP/GSNAP
*   BLAT
*   novoalign
*   BWA
*   MAQ
*   SOAP
*   **STAR**
*   **HISAT2**
*   **Kallisto**
*   Salmon

We used HISAT2 and HTSeq-count to map the sequence data from 2 samples to the reference genome in a splice-aware manner and used the Ensembl Reference Transcriptome to quantitate read levels across transcripts. HISAT2 shatters reads into fragments and to map them to the genome and splice sites.

* * *

Each run could require > 4Gb of RAM. Total run-time can be anything from 1 hour to 24 hours per sample.

You can download the Latest Human Genome build (90) from Ensembl. This is available [here](http://www.ensembl.org/info/data/ftp/index.html).

<br>
To build a new index you will need to use the hisat2-build command

We will only build an index for Chromosome 22 instead of the whole Human Genome.

```
hisat2-build Homo_sapiens.GRCh37.75.dna.chromosome.22.fa Homo_sapiens.GRCh37.dna.22
```

We also need to assemble a list of known spice sites for HISAT2, we use a utility script called _hisat2_extract_splice_sites.py_ to do this.

```
hisat2_extract_splice_sites.py Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf > known_splice_sites.txt
```

This will produce the following output file containing coordinates for all splice sites on Chr22.
```
22	16062315	16062810	+
22	16100647	16101275	-
22	16101473	16118791	-
22	16118912	16124741	-
22	16151820	16162396	-
22	16157341	16164481	+
22	16159388	16162396	-
22	16162387	16164481	+
22	16162486	16186810	-
22	16164568	16171951	+
22	16171606	16171951	+
```

The two samples (Lanes 1-4, 2 replicates) are represented by single-end sequencing files.

*   mir210_lane1.fq.gz
*   mir210_lane2.fq.gz
*   control_lane1.fq.gz
*   control_lane2.fq.gz

We will launch HISAT2 on each of the files (HISAT2 can also process sets of paired end files). By default HISAT2 assumes the sequencing was unstranded. We provide the splice site information and the HISAT2 Index. The -p 4 option asks for four processors per run, for bigger machines you can increase this for faster runs when aligning more reads.

```
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh37.dna.22 -U mir210_lane1.fq.gz | samtools view -bS - > mir210_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh37.dna.22 -U mir210_lane2.fq.gz | samtools view -bS - > mir210_lane2.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh37.dna.22 -U control_lane1.fq.gz | samtools view -bS - > control_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh37.dna.22 -U control_lane2.fq.gz | samtools view -bS - > control_lane2.bam
```

```
[samopen] SAM header is present: 194 sequences.
3243206 reads; of these:
  3243206 (100.00%) were unpaired; of these:
    3112883 (95.98%) aligned 0 times
    104048 (3.21%) aligned exactly 1 time
    26275 (0.81%) aligned >1 times
4.02% overall alignment rate
```
Only 4% aligns because we are only using chr22. You would hope/expect more to align to the whole transcriptome.

Now we can run HTSeq-Count on the results to produce count tables

**[http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)**

If the data is strand-specific you should use **-s yes**. In this case it is not strand-specific. It is also worth noting that if you are processing paired-end data you may need to sort the bam first before running _htseq-count_. This can be done using "_samtools sort -n BAM_FILE OUT_BASENAME_ ".

```
htseq-count -f bam mir210_lane1.bam Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > mir210_lane1.bam.counts
htseq-count -f bam mir210_lane2.bam Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > mir210_lane2.bam.counts
htseq-count -f bam control_lane1.bam  Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > control_lane1.bam.counts
htseq-count -f bam control_lane2.bam  Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > control_lane2.bam.counts
```

```
58105 GFF lines processed.
6160 SAM alignments  processed.
```

The Final counts tables can be loaded directly into DESeq2 in R/BioConductor. The counts are in the files with **.counts** extensions.

Lets take a look:
```
ls  *.counts
```

```
less mir210_lane1.bam.counts
```

Typically we would now switch to R/BioConductor to load these count files directly into DESeq2 for a negative binomial analysis of count statistics. This will happen in the next practical.
