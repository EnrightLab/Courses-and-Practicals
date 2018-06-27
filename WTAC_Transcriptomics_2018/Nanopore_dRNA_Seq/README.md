# Nanopore Direct RNA sequencing analysis

Adrien Leger, EMBL-EBI

26th June 2018



## Pros and cons of dRNA

#### Cons

- Lower yield than cDNA seq => 500Mb to 1.5Gb per flowcell 
- Slighly higher error ~ 11%
- Requires a massive quantity of polyA+ RNA (or of a specific target)



#### Pros

- Library Preparation much simpler, fewer steps => less biaised

- Longer reads than cDNA

  ![](pictures/Slide17.png) 

- Better exons connectivity

![](pictures/exon_align.png)

* PolyAs can be detected and measured (hopefully)

  ![](pictures/Slide09.png)

* (Some) RNA modifications can modify the signal

  ![](/home/aleg/Drive/EBI/Teaching/Courses-and-Practicals/WTAC_Transcriptomics_2018/Nanopore_dRNA_Seq/pictures/Slide34.png)

   



## ONT Fast5 file format

MINKnow generates files containing the raw intensity signal in [HDF5 format](https://support.hdfgroup.org/HDF5/). Each read is contained in a single file.

![](pictures/HDF5.jpeg)



Files can be explored using **HDFview**

**Fast5 containing raw data only**

![](pictures/fast5_pre.png)

**Fast5 after basecalling with Albacore (or live basecalling with MINKnow)**

![](pictures/fast5_post.png)



**Example of raw signal**

![](pictures/Raw1.png)

![](pictures/Raw2.png)

![](pictures/Raw3.png)



## Useful tools for dRNA-Seq ONT analysis

### Basecalling

* **Albacore** (available to ONT community)
* Guppy = GPU accelerated version of Albacore (available to developer contract with ONT)
* [Chiron](https://github.com/haotianteng/chiron) = Community alternative

See basecaller comparison => https://github.com/rrwick/Basecalling-comparison

### Quality control

* [NanoPack](https://github.com/wdecoster/nanopack) = Suite of tools to QC and process raw ONT data
* [pycoQC](https://github.com/a-slide/pycoQC) = ONT data QC for Jupyter Notebook

### RNA alignment

* [**Minimap2** ](https://github.com/lh3/minimap2)
* [STAR](https://github.com/alexdobin/STAR)
* Transcriptome based aligner => Kallisto or Salmon

### Read Polishing

* [Nanopolish event align](https://nanopolish.readthedocs.io/en/latest/)
* [Tombo resquiggle](https://nanoporetech.github.io/tombo/)

### DNA/RNA modification detection

* [Tombo detect_modifications](https://nanoporetech.github.io/tombo/)
* [Nanopolish call-methylation](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html)

Exhaustive list of tools => https://docs.google.com/spreadsheets/d/15LWXg0mUeNOHVthl8JRX-FzJ9w8jrWogS4YhDcxyAfI/edit#gid=0
