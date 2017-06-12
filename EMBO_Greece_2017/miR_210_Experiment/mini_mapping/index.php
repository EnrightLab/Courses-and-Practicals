<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<!-- saved from url=(0056)http://www.sanger.ac.uk/Teams/Team101/wtac/lesson02.html -->
<?php include ("../header_int.php"); ?>

<H1>Basic Mapping Practical</H1>
<HR>

<H4>Raw Data</H4>
<p>You can download the raw data for this practical <a href="data">here</a></p>

<H4>Initial Setup</H4>
<p>We have made a simple test with only Chromosome 22 and 400k reads on your machines</p>
<p>We will open a terminal and change to the practical folder on your desktop</p>
<pre class="com">
cd ~/Desktop/mapping_practical
</pre>

<H4>Commonly used assembly and transcriptome analysis tools</H4>
<TABLE BORDER="1">
<TR><TD BGCOLOR="lightblue">Method</TD><TD BGCOLOR="lightblue">Known Transcriptome</TD><TD BGCOLOR="lightblue">de novo Assembly</TD><TD BGCOLOR="lightblue">Ab initio Assembly</TD></TR>
<TR><TD>Cufflinks</TD><TD BGCOLOR="pink">Yes</TD><TD BGCOLOR="pink">Yes</TD><TD>No</TD></TR>
<TR><TD>Trinity</TD><TD>No</TD><TD>No</TD><TD BGCOLOR="pink">Yes</TD></TR>
<TR><TD>Scripture</TD><TD BGCOLOR="pink">Yes</TD><TD BGCOLOR="pink">Yes</TD><TD>No</TD></TR>
<TR><TD>TransAbyss</TD><TD>No</TD><TD>No</TD><TD BGCOLOR="pink">Yes</TD></TR>
<TR><TD>Velvet<TD><TD>No</TD><TD BGCOLOR="pink">Yes</TD></TR>
<TR><TD>Rum</TD><TD BGCOLOR="pink">Yes</TD><TD BGCOLOR="pink">Yes</TD><TD>No</TD></TR>
</TABLE>

<p>Many of these tools require you to perform your own mapping using a mapping algorithm</p>

<H4>Commonly used mapping tools</H4>
<UL>
<LI>BowTie
<LI>GMAP/GSNAP
<LI>BLAT
<LI>novoalign
<LI>BWA
<LI>MAQ
<LI>SOAP
<LI>Star
<LI>HiSat
<LI>Callisto
</UL>

<p>We used HiSat2 and HTSeq-count to map the sequence data from 2 samples to the reference genome in a splice-aware manner and used the Ensembl Reference Transcriptome to quantitate read levels across transcripts. 
HiSat2 shatters reads into fragments and to map them to the genome and splice sites.</p>

<HR>
<p>Each run could require > 4Gb of RAM. Total run-time can be anything from 1 hour to 24 hours per sample. </P>
</p>
 
<p>You can download the Latest Human Genome build (84) from Ensembl.
This is available <a href="http://www.ensembl.org/info/data/ftp/index.html">here</a>.</p>


<p>To actually build a new bowtie index you can use the bowtie2-build command</P>
<p>We will actually only build Chromosome 22 instead of the whole Human Genome.</p>


<pre class="com">
hisat2-build Homo_sapiens.GRCh37.75.dna.chromosome.22.fa Homo_sapiens.GRCh38.dna.22
</pre>

We also need to assemble a list of known spice sites for HiSat2, we use a utility script
called <i>hisat2_extract_splice_sites.py</i> to do this.

<pre class="com">
./hisat2_extract_splice_sites.py Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf > known_splice_sites.txt
</pre>

<p>The two samples (Lanes 1-4, 2 replicates) are represented by single end sequencing files.</p>

<UL>
<LI>mir210_lane1.fq.gz
<LI>mir210_lane2.fq.gz
<LI>control_lane1.fq.gz
<LI>control_lane2.fq.gz
</UL>

<p>We will launch HiSat2 on each of the files (HiSat can also process sets of paired end files). We also need to let TopHat know the type of sequencing (unstranded).
We provide the splice site information and the HiSat Index.
The -p 4 option asks for four processors per run, for bigger machines you can increase this for faster runs when aligning more reads.</P>

<pre class="com">
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U mir210_lane1.fq.gz | samtools view -bS - > mir210_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U mir210_lane2.fq.gz | samtools view -bS - > mir210_lane2.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U control_lane1.fq.gz | samtools view -bS - > control_lane1.bam
hisat2 --known-splicesite-infile known_splice_sites.txt -p 4 -x Homo_sapiens.GRCh38.dna.22 -U control_lane2.fq.gz | samtools view -bS - > control_lane2.bam
</pre>

<pre class="res">
[samopen] SAM header is present: 194 sequences.
3260345 reads; of these:
  3260345 (100.00%) were unpaired; of these:
    439505 (13.48%) aligned 0 times
    2302177 (70.61%) aligned exactly 1 time
    518663 (15.91%) aligned >1 times
86.52% overall alignment rate
</pre>

<p>Now we can run HTSeq-Count on the tophat results to produce count tables</p>

<b><a href="http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html">http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html</a></b>

<p>If the data is strand-specific you should use <b>-s yes</b>. In this case it is not strand-specific. It is also worth noting that if you are processing paired-end data 
you may need to sort the bam first before running <i>htseq-count</i>. This can be done using "<i>samtools sort -n  BAM_FILE  OUT_BASENAME </i>".</p>

<pre class="com">
htseq-count -f bam mir210_lane1.bam Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > mir210_lane1.bam.counts
htseq-count -f bam mir210_lane2.bam Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > mir210_lane2.bam.counts
htseq-count -f bam control_lane1.bam  Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > control_lane1.bam.counts
htseq-count -f bam control_lane2.bam  Homo_sapiens.GRCh37.75.dna.chromosome.22.gtf --stranded=no > control_lane2.bam.counts
</pre>

<pre class="res">
58105 GFF lines processed.
6160 SAM alignments  processed.
</pre>

<p>The Final counts tables can be loaded directly into DESeq2 in R/BioConductor
the counts are in the files with <b>.counts</b> extensions</P>

</html>
