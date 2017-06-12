Mini-mapping Practical
================
Anton Enright & Dimitrios Vitsios
'12 June, 2017'

Basic Mapping Practical
-------

## Raw Data

You can download the raw data for this practical [here](data)

## Initial Setup

We have made a simple test with only Chromosome 22 and 400k reads on your machines

We will open a terminal and change to the practical folder on your desktop

```
cd ~/Desktop/mapping_practical
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

*   BowTie
*   GMAP/GSNAP
*   BLAT
*   novoalign
*   BWA
*   MAQ
*   SOAP
*   Star
*   HiSat
*   Callisto

We used HiSat2 and HTSeq-count to map the sequence data from 2 samples to the reference genome in a splice-aware manner and used the Ensembl Reference Transcriptome to quantitate read levels across transcripts. HiSat2 shatters reads into fragments and to map them to the genome and splice sites.

* * *

Each run could require > 4Gb of RAM. Total run-time can be anything from 1 hour to 24 hours per sample.

You can download the Latest Human Genome build (84) from Ensembl. This is available [here](http://www.ensembl.org/info/data/ftp/index.html).

To actually build a new bowtie index you can use the bowtie2-build command

We will actually only build Chromosome 22 instead of the whole Human Genome.

```
hisat2-build Homo_sapiens.GRCh37.75.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.22
```

We also need to assemble a list of known spice sites for HiSat2, we use a utility script called _hisat2_extract_splice_sites.py_ to do this.

```
./hisat2_extract_splice_sites.py Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf > known_splice_sites.txt
```

The two samples (Lanes 1-4, 2 replicates) are represented by single end sequencing files.

*   mir210_lane1.fq.gz
*   mir210_lane2.fq.gz
*   control_lane1.fq.gz
*   control_lane2.fq.gz

We will launch HiSat2 on each of the files (HiSat can also process sets of paired end files). We also need to let TopHat know the type of sequencing (unstranded). We provide the splice site information and the HiSat Index. The -p 4 option asks for four processors per run, for bigger machines you can increase this for faster runs when aligning more reads.

```
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U mir210_lane1.fq.gz | samtools view -bS - > mir210_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U mir210_lane2.fq.gz | samtools view -bS - > mir210_lane2.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U control_lane1.fq.gz | samtools view -bS - > control_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U control_lane2.fq.gz | samtools view -bS - > control_lane2.bam
```

```
[samopen] SAM header is present: 194 sequences.
3260345 reads; of these:
  3260345 (100.00%) were unpaired; of these:
    439505 (13.48%) aligned 0 times
    2302177 (70.61%) aligned exactly 1 time
    518663 (15.91%) aligned >1 times
86.52% overall alignment rate
```

Now we can run HTSeq-Count on the tophat results to produce count tables

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

The Final counts tables can be loaded directly into DESeq2 in R/BioConductor the counts are in the files with **.counts** extensions
